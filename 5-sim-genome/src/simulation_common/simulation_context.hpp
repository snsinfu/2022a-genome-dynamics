#pragma once

#include <md.hpp>


// Structure: simulation_context
//
// Aggregate of simulation-global variables and statistics. Used to pass the
// said values to `simulation_store`.
//
struct simulation_context
{
    // Global variables
    md::scalar time          = 0;
    md::vector wall_semiaxes = {0, 0, 0};
    md::scalar bead_scale    = 1;
    md::scalar bond_scale    = 1;

    // Statistics
    md::scalar mean_energy = 0;
    md::scalar wall_energy = 0;
};
