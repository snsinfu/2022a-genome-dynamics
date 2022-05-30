#pragma once

#include <string>

#include <h5.hpp>
#include <md.hpp>

#include "simulation_config.hpp"
#include "simulation_data.hpp"


class simulation_store
{
public:
    explicit
    simulation_store(std::string const& filename);

    void save_config(simulation_config const& config);
    void save_beads(md::array_view<bead_data const> beads);
    void save_chains(md::array_view<chain_data const> chains);
    void save_snapshot(md::step step, md::array_view<md::point const> positions);

private:
    h5::file _file;
};
