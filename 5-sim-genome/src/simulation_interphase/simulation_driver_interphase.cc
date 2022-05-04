#include <cmath>

#include <md.hpp>

#include "simulation_driver.hpp"


void simulation_driver::run_simulation()
{
    _store.set_phase("interphase");

    auto callback = [=](md::step step) {
        _context.time = step * _config.interphase_timestep;

        // Calculating energy is expensive. So update stats only when needed.
        auto const with_logging = step % _config.interphase_logging_interval == 0;
        auto const with_sampling = step % _config.interphase_sampling_interval == 0;
        auto const sample_frame = step / _config.interphase_sampling_interval;

        if (with_logging || with_sampling) {
            _context.mean_energy = _system.compute_energy() / _system.particle_count();
        }

        if (with_logging) {
            print_progress("inter", step);
        }

        if (with_sampling) {
            _store.save_positions(step, _system.view_positions());
            _store.save_context(step, _context);
        }

        if (step % _config.contactmap_update_interval == 0) {
            _contact_map.update(_system.view_positions());
        }

        if (with_sampling && sample_frame % _config.contactmap_thinning_rate == 0) {
            _store.save_contacts(step, _contact_map.accumulate());
            _contact_map.clear();
        }

        update_bead_scale();
        update_wall_semiaxes();
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.interphase_temperature,
        .spacestep   = _config.interphase_spacestep,
        .timestep    = _config.interphase_timestep,
        .steps       = _config.interphase_steps,
        .seed        = _random(),
        .callback    = callback
    });
}


void simulation_driver::update_bead_scale()
{
    auto const t_bead = _context.time / _config.bead_scale_tau;
    auto const t_bond = _context.time / _config.bond_scale_tau;
    _context.bead_scale = 1 - (1 - _config.bead_scale_init) * std::exp(-t_bead);
    _context.bond_scale = 1 - (1 - _config.bond_scale_init) * std::exp(-t_bond);

    _contact_map.set_contact_distance(_config.contactmap_distance * _context.bead_scale);
}


void simulation_driver::update_wall_semiaxes()
{
    auto& semiaxes = _context.wall_semiaxes;

    md::vector net_force;
    net_force += _compute_packing_reaction();
    net_force -= _config.wall_semiaxes_spring.hadamard(semiaxes);

    // Simulate ad-hoc overdamped motion of the wall.
    semiaxes += _config.interphase_timestep * _config.wall_mobility * net_force;
}
