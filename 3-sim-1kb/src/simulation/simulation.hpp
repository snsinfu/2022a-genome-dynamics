#pragma once

#include <memory>
#include <random>

#include <md.hpp>

#include "config.hpp"
#include "glues.hpp"
#include "loops.hpp"
#include "store.hpp"
#include "topology.hpp"


class simulation
{
public:
    explicit simulation(simulation_config const& config);
    void run();

private:
    void setup_particles();
    void setup_loops();
    void setup_glues();
    void setup_forcefield_repulsion();
    void setup_forcefield_connectivity();
    void setup_forcefield_bending();
    void setup_forcefield_loop();
    void setup_forcefield_glue();
    void initialize_particles();
    void initialize_loops();
    void run_simulation();
    void step_loops(md::step step);
    void step_glues(md::step step);
    void show_progress(md::step step);
    void save_sample();

private:
    simulation_config               _config;
    md::system                      _system;
    std::vector<chain_assignment>   _chains;
    std::shared_ptr<loop_simulator> _loops;
    std::shared_ptr<glue_simulator> _glues;
    std::mt19937_64                 _random;
    simulation_store                _store;
};
