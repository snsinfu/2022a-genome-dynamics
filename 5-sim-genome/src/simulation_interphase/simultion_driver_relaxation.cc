#include <algorithm>

#include <md.hpp>

#include "simulation_driver.hpp"


void simulation_driver::run_relaxation()
{
    _store.set_phase("relaxation");
    {
        auto const init_positions = _store.load_positions(0);
        auto positions = _system.view_positions();
        std::copy(init_positions.begin(), init_positions.end(), positions.begin());
    }

    auto callback = [=](md::step step) {
        // Calculating energy is expensive. So update stats only when needed.
        auto const with_logging = step % _config.relaxation_logging_interval == 0;
        auto const with_sampling = step % _config.relaxation_sampling_interval == 0;

        if (with_logging || with_sampling) {
            _context.mean_energy = _system.compute_energy() / _system.particle_count();
        }

        if (with_logging) {
            print_progress("relax", step);
        }

        if (with_sampling) {
            _store.save_positions(step, _system.view_positions());
            _store.save_context(step, _context);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(_system, {
        .temperature = _config.relaxation_temperature,
        .spacestep   = _config.relaxation_spacestep,
        .timestep    = _config.relaxation_timestep,
        .steps       = _config.relaxation_steps,
        .seed        = _random(),
        .callback    = callback
    });
}
