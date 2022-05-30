#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>

#include <md.hpp>

#include "simulation_config.hpp"
#include "simulation_data.hpp"
#include "simulation_driver.hpp"
#include "simulation_store.hpp"
#include "walltime.hpp"


namespace
{
    struct particle_data
    {
        md::scalar a_factor = 0;
        md::scalar b_factor = 0;
    };

    md::attribute_key<particle_data> particle_data_attribute;
}


simulation_driver::simulation_driver(simulation_config const& config)
    : _config{config}
    , _store{config.output}
    , _random{config.seed}
{
    _store.save_config(_config);
    setup_particles();
    setup_forcefield();
}


void simulation_driver::setup_particles()
{
    auto const beads = load_beads_data(_config.beads_filename);

    std::string cur_chain = "";
    md::index cur_start = 0;
    md::index cur_end = 0;

    auto finish_chain = [&] {
        if (cur_start == cur_end) {
            return;
        }
        _chains.push_back({
            .start = cur_start,
            .end   = cur_end,
        });
    };

    _system.require(particle_data_attribute);

    for (auto const& bead : beads) {
        auto part = _system.add_particle({
            .mobility = _config.mobility,
        });
        part.view(particle_data_attribute) = {
            .a_factor = bead.a_factor,
            .b_factor = bead.b_factor,
        };
        cur_end = part.index;

        if (bead.chain != cur_chain) {
            finish_chain();
            cur_chain = bead.chain;
            cur_start = cur_end;
        }
    }

    // Make cur_end past the end of the last chain.
    cur_end++;
    finish_chain();

    // Save beads and topology to trajectory.
    _store.save_beads(beads);
    _store.save_chains(_chains);
}


void simulation_driver::setup_forcefield()
{
    setup_forcefield_repulsions();
    setup_forcefield_bonds();
    setup_forcefield_outer_wall();
    setup_forcefield_inner_wall();
}


void simulation_driver::setup_forcefield_repulsions()
{
    // Short-range mixed repulsion.

    md::softcore_potential<2> const a_potential {
        .energy   = _config.a_core_repulsion,
        .diameter = _config.a_core_diameter,
    };

    md::softcore_potential<8> const b_potential {
        .energy   = _config.b_core_repulsion,
        .diameter = _config.b_core_diameter,
    };

    _system.add_forcefield(
        md::make_neighbor_pairwise_forcefield(
            [=](md::index i, md::index j) {
                auto const data = _system.view(particle_data_attribute);
                auto const a_factor = (data[i].a_factor + data[j].a_factor) / 2;
                auto const b_factor = (data[i].b_factor + data[j].b_factor) / 2;
                return a_factor * a_potential + b_factor * b_potential;
            }
        )
        .set_neighbor_distance(
            std::max(a_potential.diameter, b_potential.diameter)
        )
    );
}


void simulation_driver::setup_forcefield_bonds()
{
    auto bonds = _system.add_forcefield(
        md::make_bonded_pairwise_forcefield(
            md::harmonic_potential {
                .spring_constant = _config.bond_spring,
            }
        )
    );

    for (auto const& chain : _chains) {
        bonds->add_bonded_range(chain.start, chain.end);
    }
}


void simulation_driver::setup_forcefield_outer_wall()
{
    md::sphere const outer_wall {
        .radius = _config.outer_wall_radius
    };
    md::scalar const wall_a_factor = 0;
    md::scalar const wall_b_factor = 1;

    md::softcore_potential<2> const a_potential {
        .energy   = _config.a_core_repulsion,
        .diameter = _config.a_core_diameter / 2,
    };

    md::softcore_potential<8> const b_potential {
        .energy   = _config.b_core_repulsion,
        .diameter = _config.b_core_diameter / 2,
    };

    _system.add_forcefield(
        md::make_sphere_inward_forcefield(
            [=](md::index i) {
                auto const data = _system.view(particle_data_attribute);
                auto const a_factor = (data[i].a_factor + wall_a_factor) / 2;
                auto const b_factor = (data[i].b_factor + wall_b_factor) / 2;
                return _config.outer_wall_multiplier * (
                    a_factor * a_potential + b_factor * b_potential
                );
            }
        )
        .set_sphere(outer_wall)
    );

    _system.add_forcefield(
        md::make_sphere_outward_forcefield(
            md::harmonic_potential {
                .spring_constant = _config.outer_wall_spring
            }
        )
        .set_sphere(outer_wall)
    );
}


void simulation_driver::setup_forcefield_inner_wall()
{
    if (_config.inner_wall_radius < 1e-6) {
        return;
    }

    md::sphere const inner_wall {
        .radius = _config.inner_wall_radius
    };
    md::scalar const wall_a_factor = 0;
    md::scalar const wall_b_factor = 1;

    md::softcore_potential<2> const a_potential {
        .energy   = _config.a_core_repulsion,
        .diameter = _config.a_core_diameter / 2,
    };

    md::softcore_potential<8> const b_potential {
        .energy   = _config.b_core_repulsion,
        .diameter = _config.b_core_diameter / 2,
    };

    _system.add_forcefield(
        md::make_sphere_inward_forcefield(
            md::harmonic_potential {
                .spring_constant = _config.inner_wall_spring
            }
        )
        .set_sphere(inner_wall)
    );

    _system.add_forcefield(
        md::make_sphere_outward_forcefield(
            [=](md::index i) {
                auto const data = _system.view(particle_data_attribute);
                auto const a_factor = (data[i].a_factor + wall_a_factor) / 2;
                auto const b_factor = (data[i].b_factor + wall_b_factor) / 2;
                return _config.inner_wall_multiplier * (
                    a_factor * a_potential + b_factor * b_potential
                );
            }
        )
        .set_sphere(inner_wall)
    );
}


void simulation_driver::run()
{
    run_initialization();
    run_sampling();
}


void simulation_driver::run_initialization()
{
    auto positions = _system.view_positions();

    for (auto const& chain : _chains) {
        std::uniform_real_distribution<md::scalar> center_coord {
            -_config.outer_wall_radius, _config.outer_wall_radius
        };
        md::point center;
        do {
            center = {
                center_coord(_random), center_coord(_random), center_coord(_random),
            };
        } while (center.distance({0, 0, 0}) < _config.inner_wall_radius);

        std::normal_distribution<md::scalar> normal;
        md::vector const direction = md::normalize({
                normal(_random), normal(_random), normal(_random),
        });

        md::vector delta;
        md::point pos;
        for (md::index i = chain.start; i < chain.end; i++) {
            positions[i] = pos;
            delta += pos - center;
            pos += _config.init_bond_length * direction;
        }
        delta /= md::scalar(chain.end - chain.start);

        // Correct the centroid
        for (md::index i = chain.start; i < chain.end; i++) {
            positions[i] -= delta;
        }
    }
}


void simulation_driver::run_sampling()
{
    std::clog << "[sim] sampling...\n";

    auto const log_progress = [&](md::step step) {
        auto const energy = _system.compute_energy() / md::scalar(_system.particle_count());
        std::clog
            << "[sim] "
            << walltime::now()
            << '\t'
            << step
            << '\t'
            << "E: "
            << energy
            << '\n';
    };

    auto const callback = [&](md::step step) {
        if (step % _config.logging_interval == 0) {
            log_progress(step);
        }
        if (step % _config.sampling_interval == 0) {
            _store.save_snapshot(step, _system.view_positions());
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.temperature,
        .timestep    = _config.timestep,
        .steps       = _config.steps,
        .callback    = callback,
    });
}
