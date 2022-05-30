#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include <md.hpp>

#include "../simulation_common/particle_data.hpp"
#include "../simulation_common/simulation_config.hpp"
#include "../simulation_common/simulation_store.hpp"

#include "simulation_driver.hpp"


namespace
{
    template<typename T>
    std::shared_ptr<T> copy_shared(T const& obj)
    {
        return std::make_shared<T>(obj);
    }
}


simulation_driver::simulation_driver(simulation_store& store)
    : _store{store}
    , _config{store.load_config()}
    , _random{_config.spindle_seed}
{
    setup();
}


void simulation_driver::setup()
{
    setup_chains();
    setup_particles();
    setup_forcefield();
}


void simulation_driver::setup_chains()
{
    auto const chroms = _store.load_chromosomes();
    auto const coarse = _config.init_coarse_graining;

    md::index start = 0;

    for (auto const& chrom : chroms) {
        auto const size = chrom.end - chrom.start;
        auto const cen = (chrom.centromere_start + chrom.centromere_end) / 2;
        auto const coarse_size = (size + coarse - 1) / coarse;
        auto const coarse_cen = (cen - chrom.start) / coarse;

        _chains.push_back({
            .chromosome = chrom,
            .start      = start,
            .end        = start + coarse_size,
            .centromere = start + coarse_cen
        });
        start += coarse_size;
    }
}


void simulation_driver::setup_particles()
{
    for (auto const& chain : _chains) {
        for (md::index i = chain.start; i < chain.end; i++) {
            _system.add_particle({
                .mobility = _config.init_mobility
            });
        }
    }
}


void simulation_driver::setup_forcefield()
{
    setup_repulsion_forcefield();
    setup_connectivity_forcefield();
    setup_spindle_forcefield();
    setup_packing_forcefield();
}


void simulation_driver::setup_repulsion_forcefield()
{
    // General repulsion for avoiding chain crossings.

    _system.add_forcefield(
        md::make_neighbor_pairwise_forcefield(
            md::softcore_potential<2, 3> {
                .energy   = _config.init_bead_repulsion,
                .diameter = _config.init_bead_diameter
            }
        )
        .set_neighbor_distance(_config.init_bead_diameter)
    );
}


void simulation_driver::setup_connectivity_forcefield()
{
    // Spring bonds and bending cost.

    auto bonds = _system.add_forcefield(
        md::make_bonded_pairwise_forcefield(
            md::semispring_potential {
                .spring_constant      = _config.init_bond_spring,
                .equilibrium_distance = _config.init_bond_length
            }
        )
    );

    auto bends = _system.add_forcefield(
        md::make_bonded_triplewise_forcefield(
            md::cosine_bending_potential {
                .bending_energy = _config.init_bend_energy
            }
        )
    );

    for (auto const& chain : _chains) {
        bonds->add_bonded_range(chain.start, chain.end);
        bends->add_bonded_range(chain.start, chain.end);
    }
}


void simulation_driver::setup_spindle_forcefield()
{
    // Spindle core attracts centromeres.

    std::vector<md::index> centromeres;
    for (auto const& chain : _chains) {
        assert(chain.centromere > chain.start);
        assert(chain.centromere + 1 < chain.end);
        centromeres.push_back(chain.centromere - 1);
        centromeres.push_back(chain.centromere);
        centromeres.push_back(chain.centromere + 1);
    }

    _spindle_forcefield = copy_shared(
        md::make_point_source_forcefield(
            md::harmonic_potential {
                .spring_constant = _config.init_spindle_spring
            }
        )
        .set_point_source(_config.init_spindle_point)
        .set_point_source_targets(centromeres)
    );
}


void simulation_driver::setup_packing_forcefield()
{
    // Weak harmonic well potential prevents open diffusion.

    _packing_forcefield = copy_shared(
        md::make_point_source_forcefield(
            md::semispring_potential {
                .spring_constant      = _config.init_packing_spring,
                .equilibrium_distance = _config.init_packing_radius
            }
        )
        .set_point_source(_config.init_spindle_point)
    );
}


void simulation_driver::run()
{
    run_initialization();
    run_spindle_phase();
    run_packing_phase();
}


void simulation_driver::run_initialization()
{
    // Initialize chains as randomly-directed rods.
    auto positions = _system.view_positions();

    for (auto const& chain : _chains) {
        std::normal_distribution<md::scalar> normal;

        auto const centroid = _config.init_start_point + _config.init_start_stddev * md::vector {
            normal(_random), normal(_random), normal(_random)
        };
        auto const step = _config.init_bond_length * md::normalize(md::vector {
            normal(_random), normal(_random), normal(_random)
        });

        auto pos = centroid - step * (chain.end - chain.start) / 2;
        for (md::index i = chain.start; i < chain.end; i++) {
            positions[i] = pos;
            pos += step;
        }
    }
}


void simulation_driver::run_spindle_phase()
{
    _store.set_phase("spindle");
    save_chains();

    _system.add_forcefield(_spindle_forcefield);

    auto const callback = [&](md::step step) {
        if (step % _config.init_sampling_interval == 0) {
            _store.save_positions(step, _system.view_positions());
        }

        if (step % _config.init_logging_interval == 0) {
            print_progress("spindle", step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.init_temperature,
        .spacestep   = _config.init_spacestep,
        .timestep    = _config.init_timestep,
        .steps       = _config.init_spindle_steps,
        .seed        = _random(),
        .callback    = callback
    });
}


void simulation_driver::run_packing_phase()
{
    _store.set_phase("packing");
    save_chains();

    _system.add_forcefield(_packing_forcefield);

    auto const callback = [&](md::step step) {
        if (step % _config.init_sampling_interval == 0) {
            _store.save_positions(step, _system.view_positions());
        }

        if (step % _config.init_logging_interval == 0) {
            print_progress("packing", step);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.init_temperature,
        .spacestep   = _config.init_spacestep,
        .timestep    = _config.init_timestep,
        .steps       = _config.init_packing_steps,
        .seed        = _random(),
        .callback    = callback
    });
}


void simulation_driver::print_progress(std::string const& phase, md::step step)
{
    auto const wallclock_time = std::time(nullptr);

    std::clog
        << "[" + phase + "] "
        << std::put_time(std::localtime(&wallclock_time), "%F %T")
        << '\t'
        << step
        << '\t'
        << "E: "
        << _system.compute_energy() / _system.particle_count()
        << '\n';
}


void simulation_driver::save_chains()
{
    std::vector<chromosome_range> chroms;

    for (auto const& chain : _chains) {
        chroms.push_back({
            .name  = chain.chromosome.name,
            .start = chain.start,
            .end   = chain.end
        });
    }

    _store.save_chromosomes(chroms);
}
