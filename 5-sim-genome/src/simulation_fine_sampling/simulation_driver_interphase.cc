#include <cmath>

#include <md.hpp>

#include "simulation_driver.hpp"


void simulation_driver::run_simulation()
{
    _store.set_phase("fine_sampling");

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
            print_progress("fine", step);
        }

        if (with_sampling) {
            _store.save_positions(step, _system.view_positions());
            _store.save_context(step, _context);
        }

        update_wall_semiaxes();
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.interphase_temperature,
        .timestep    = _config.interphase_timestep,
        .spacestep   = _config.interphase_spacestep,
        .steps       = _config.interphase_steps,
        .seed        = _random(),
        .callback    = callback
    });
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
