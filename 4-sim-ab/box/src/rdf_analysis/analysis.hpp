#pragma once

#include <string>

#include <md.hpp>


struct analysis_config
{
    std::string filename;

    std::string particle_selection;
    md::step step_start = 0;
    md::step step_end = 0;
    md::scalar bin_width = 0.1;
    md::scalar max_distance = 1;
};

void run_analysis(analysis_config const& config);
