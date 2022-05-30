#pragma once

#include <md.hpp>

#include "simulation_config.hpp"
#include "simulation_data.hpp"
#include "simulation_store.hpp"


// Manages a single simulation run.
class simulation_driver
{
public:
    // Configures a simulation system with the given config parameters.
    explicit simulation_driver(simulation_config const& config);

    // Runs a single complete simulation: initialization, relaxation and
    // sampling steps. Writes snapshots to a trajectory file.
    void run();

private:
    void setup_particles();
    void setup_forcefield();
    void setup_forcefield_repulsions();
    void setup_forcefield_bonds();
    void run_initialization();
    void run_sampling();

private:
    simulation_config _config;
    simulation_store _store;
    md::system _system;
    md::random_engine _random;
    std::vector<chain_data> _chains;
};
