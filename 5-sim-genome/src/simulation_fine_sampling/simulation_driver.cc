#include <algorithm>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>

#include <md.hpp>

#include "../simulation_common/simulation_store.hpp"

#include "simulation_driver.hpp"


simulation_driver::simulation_driver(simulation_store& store)
    : _store{store}
    , _config{store.load_config()}
    , _random{_config.interphase_seed ^ 700000} // ?
{
    setup();
}


void simulation_driver::setup()
{
    setup_particles();
    setup_forcefield();
    setup_context();

    // ?
    _config.interphase_temperature = 0;
    _config.interphase_sampling_interval = 100;
    _config.interphase_steps = 1000 * 100;
    _config.interphase_timestep = 1e-5 / 100;
}


void simulation_driver::setup_context()
{
    // Load from trajectory file.

    _store.set_phase("interphase");

    auto const step = 700000;
    auto const init_positions = _store.load_positions(step);
    auto const init_context = _store.load_context(step);
    auto positions = _system.view_positions();
    std::copy(init_positions.begin(), init_positions.end(), positions.begin());

    _context = init_context;

    // ?
    _context.bead_scale = 1;
    _context.bond_scale = 1;
}


void simulation_driver::run()
{
    run_simulation();
}


void simulation_driver::print_progress(std::string phase, md::step step)
{
    auto const wallclock_time = std::time(nullptr);
    auto const effective_radius = std::cbrt(
        _context.wall_semiaxes.x *
        _context.wall_semiaxes.y *
        _context.wall_semiaxes.z
    );

    std::clog
        << "[" + phase + "] "
        << std::put_time(std::localtime(&wallclock_time), "%F %T")
        << '\t'
        << step
        << '\t'
        << "t: "
        << _context.time
        << '\t'
        << "R: "
        << effective_radius
        << '\t'
        << "E: "
        << _context.mean_energy
        << '\n';
}


void simulation_driver::save_chains()
{
    auto const chroms = _store.load_chromosomes();
    _store.save_chromosomes(chroms);
}
