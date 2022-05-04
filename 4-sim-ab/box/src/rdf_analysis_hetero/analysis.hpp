#pragma once

#include <string>

#include <md.hpp>


struct analysis_config
{
    std::string filename;
    std::string center_type;
    md::scalar bin_width = 0.1;
    md::scalar max_distance = 1;
};

void run_analysis(analysis_config const& config);
