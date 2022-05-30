#pragma once

#include <vector>

#include <md.hpp>

#include "config.hpp"


struct chain_assignment
{
    md::index    start = 0;
    md::index    end   = 0;
    chain_config config;
};


/**
 * Given simulation configuration, computes the assignments of indices into a
 * flat array of chain monomers.
 */
std::vector<chain_assignment> make_chain_assignments(simulation_config const& config);
