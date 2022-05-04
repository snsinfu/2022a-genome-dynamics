#pragma once

#include <memory>
#include <random>
#include <string>

#include <md.hpp>

#include "../simulation_common/particle_data.hpp"
#include "../simulation_common/simulation_config.hpp"
#include "../simulation_common/simulation_store.hpp"


struct coarse_chain_range
{
    chromosome_range chromosome;
    md::index        start;
    md::index        end;
    md::index        centromere;
};


class simulation_driver
{
public:
    explicit simulation_driver(simulation_store& store);
    void run();

private:
    void setup();
    void setup_chains();
    void setup_particles();
    void setup_forcefield();
    void setup_repulsion_forcefield();
    void setup_connectivity_forcefield();
    void setup_spindle_forcefield();
    void setup_packing_forcefield();

    void run_initialization();
    void run_spindle_phase();
    void run_packing_phase();

    void print_progress(std::string const& phase, md::step step);

    void save_chains();

private:
    simulation_store&               _store;
    simulation_config               _config;
    md::system                      _system;
    std::mt19937_64                 _random;
    std::vector<coarse_chain_range> _chains;
    std::shared_ptr<md::forcefield> _spindle_forcefield;
    std::shared_ptr<md::forcefield> _packing_forcefield;
};
