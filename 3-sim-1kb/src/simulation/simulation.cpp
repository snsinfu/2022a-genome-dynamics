#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include <md.hpp>

#include "glues.hpp"
#include "inits.hpp"
#include "loops.hpp"
#include "simulation.hpp"
#include "topology.hpp"
#include "forces/glue_forcefield.hpp"
#include "forces/loop_forcefield.hpp"


struct monomer_data
{
    md::scalar bending_energy = 0;
};

static md::attribute_key<monomer_data> monomer_data_attribute;


static std::mt19937_64
make_random(std::uint64_t seed)
{
    std::seed_seq seed_seq{seed};
    return std::mt19937_64{seed_seq};
}


simulation::simulation(simulation_config const& config)
    : _config{config}
    , _random{make_random(config.sampling.random_seed)}
    , _store{config.sampling.output_filename}
{
    setup_particles();
    setup_loops();
    setup_glues();
    setup_forcefield_repulsion();
    setup_forcefield_connectivity();
    setup_forcefield_bending();
    setup_forcefield_loop();
    setup_forcefield_glue();

    _store.save_metadata({.config = _config});
}


void
simulation::setup_particles()
{
    _system.require(monomer_data_attribute);

    _chains = make_chain_assignments(_config);

    for (auto const& chain : _chains) {
        for (md::index i = 0; i < chain.config.length; i++) {
            auto part = _system.add_particle();
            part.mobility = _config.chain.monomer_mobility;
            part.view(monomer_data_attribute) = {
                .bending_energy = _config.chain.bending_energy,
            };
        }

        // Block overrides
        auto data = _system.view(monomer_data_attribute);

        for (auto const& block : chain.config.blocks) {
            for (md::index i = block.start; i < block.end; i++) {
                if (block.bending_energy) {
                    data[chain.start + i].bending_energy = *block.bending_energy;
                }
            }
        }
    }
}


void
simulation::setup_loops()
{
    _loops = make_loop_simulator(_config);
}


void
simulation::setup_glues()
{
    _glues = make_glue_simulator(_config);
}


void
simulation::setup_forcefield_repulsion()
{
    md::softcore_potential<2, 3> const repulsive_potential {
        .energy   = _config.chain.repulsive_energy,
        .diameter = _config.chain.repulsive_diameter,
    };

    md::softcore_potential<8, 3> const attractive_potential {
        .energy   = _config.chain.attractive_energy * -1,
        .diameter = _config.chain.attractive_diameter,
    };

    auto const interaction = repulsive_potential + attractive_potential;
    auto const max_diameter = std::max(repulsive_potential.diameter, attractive_potential.diameter);
    auto const box_size = _config.chain.box_size;

    _system.add_forcefield(
        md::make_neighbor_pairwise_forcefield<md::periodic_box>(interaction)
        .set_unit_cell({
            .x_period = box_size,
            .y_period = box_size,
            .z_period = box_size,
        })
        .set_neighbor_distance(max_diameter)
    );
}


void
simulation::setup_forcefield_connectivity()
{
    md::spring_potential const potential = {
        .spring_constant      = _config.chain.bond_spring,
        .equilibrium_distance = _config.chain.bond_length,
    };

    auto bonds = _system.add_forcefield(
        md::make_bonded_pairwise_forcefield(potential)
    );

    for (auto const& chain : _chains) {
        bonds->add_bonded_range(chain.start, chain.end);
    }
}


void
simulation::setup_forcefield_bending()
{
    auto const potential = [=](md::system const& system, md::index, md::index i, md::index) {
        auto const data = system.view(monomer_data_attribute);
        return md::cosine_bending_potential{data[i].bending_energy};
    };

    auto bends = _system.add_forcefield(md::make_bonded_triplewise_forcefield(potential));

    for (auto const& chain : _chains) {
        bends->add_bonded_range(chain.start, chain.end);
    }
}


void
simulation::setup_forcefield_loop()
{
    _system.add_forcefield(
        loop_forcefield{_loops, _config.loop.bond_spring, _config.chain.repulsive_diameter}
    );
}


void
simulation::setup_forcefield_glue()
{
    auto const box_size = _config.chain.box_size;

    _system.add_forcefield(
        glue_forcefield{
            _glues,
            {.x_period = box_size, .y_period = box_size, .z_period = box_size},
            _config.glue.glue_energy,
            _config.glue.glue_distance
        }
    );
}


void
simulation::run()
{
    initialize_particles();
    initialize_loops();
    run_simulation();
}


void
simulation::initialize_particles()
{
    make_initializer(_config)->initialize(_system, _random);
}


void
simulation::initialize_loops()
{
    if (_config.sampling.loop_preloading) {
        _loops->preload(_random);
    }
}


void
simulation::run_simulation()
{
    auto const callback = [=](md::step step) {
        if (step % _config.sampling.logging_interval == 0) {
            show_progress(step);
        }

        if (step % _config.sampling.sampling_interval == 0) {
            save_sample();
        }

        if (step % _config.sampling.loop_update_interval == 0) {
            step_loops(step);
        }

        if (step % _config.sampling.glue_update_interval == 0) {
            step_glues(step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.sampling.temperature,
        .timestep    = _config.sampling.timestep,
        .steps       = _config.sampling.steps,
        .seed        = _random(),
        .callback    = callback,
    });
}


void
simulation::step_loops(md::step step)
{
    auto const leap_time =
        _config.sampling.timestep * md::scalar(_config.sampling.loop_update_interval);

    if (!_config.sampling.clear_loops_at || step < _config.sampling.clear_loops_at) {
        _loops->step(leap_time, _random);
    }

    if (step + 1 == _config.sampling.clear_loops_at) {
        _loops->clear();
    }
}


void
simulation::step_glues(md::step step)
{
    auto const leap_time =
        _config.sampling.timestep * md::scalar(_config.sampling.glue_update_interval);

    _glues->update(leap_time, _system.view_positions(), _random);
}


void
simulation::show_progress(md::step step)
{
    auto const particle_average = [&](auto number) {
        return md::scalar(number) / md::scalar(_system.particle_count());
    };

    auto const energy_metric = particle_average(_system.compute_energy());
    auto const loop_metric = particle_average(
        std::count_if(_loops->begin(), _loops->end(), [](auto const& loop) {
            return loop.id > 0;
        })
    );
    auto const glue_metric = particle_average(_glues->size());

    std::clog
        << step
        << '\t'
        << "E: " << energy_metric
        << '\t'
        << "L: " << loop_metric
        << '\t'
        << "G: " << glue_metric
        << '\n';
}


void
simulation::save_sample()
{
    _store.save_snapshot({
        .positions = _system.view_positions(),
        .loops     = md::array_view<loop_pair const>{_loops->begin(), _loops->size()},
    });
}
