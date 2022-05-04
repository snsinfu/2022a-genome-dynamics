#pragma once

#include <functional>
#include <random>

#include <md.hpp>

#include "../simulation_common/simulation_config.hpp"
#include "../simulation_common/simulation_context.hpp"
#include "../simulation_common/simulation_store.hpp"


class simulation_driver
{
public:
    explicit simulation_driver(simulation_store& store);
    void run();

private:
    void setup();
    void setup_particles();
    void setup_forcefield();
    void setup_repulsive_forcefield();
    void setup_connectivity_forcefield();
    void setup_nucleolus_forcefield();
    void setup_membrane_forcefield();
    void setup_context();

    void run_simulation();

    void print_progress(std::string phase, md::step step);

    void update_wall_semiaxes();

    void save_chains();

private:
    simulation_store&  _store;
    simulation_config  _config;
    simulation_context _context;
    md::system         _system;
    std::mt19937_64    _random;

    std::function<md::vector()> _compute_packing_reaction;
};
