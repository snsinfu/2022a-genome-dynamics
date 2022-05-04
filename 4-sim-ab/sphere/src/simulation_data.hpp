#pragma once

#include <string>
#include <vector>

#include <md.hpp>


struct bead_data
{
    std::string chain;
    md::scalar a_factor;
    md::scalar b_factor;
};


struct chain_data
{
    std::string name;
    md::index start;
    md::index end;
};


struct context_data
{
    md::step step;
    md::scalar mean_energy;
};


std::vector<bead_data> load_beads_data(std::string const& filename);
