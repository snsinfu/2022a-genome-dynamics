#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <highfive/H5File.hpp>
#include <md.hpp>

#include "particle_data.hpp"
#include "simulation_config.hpp"
#include "simulation_context.hpp"


namespace H5 = HighFive;


struct chromosome_range
{
    std::string name;
    md::index start;
    md::index end;
    md::index centromere_start;
    md::index centromere_end;
};


struct nucleolus_bond
{
    md::index nor_index;
    md::index nuc_index;
};


struct index_range
{
    md::index begin = 0;
    md::index end = 0;
};


class simulation_store
{
public:
    // Constructor takes the filename of the HDF5 file to operate on and opens
    // it in read-write mode.
    explicit simulation_store(std::string const& filename);

    // Metadata
    simulation_config             load_config();
    std::vector<chromosome_range> load_chromosomes();
    std::vector<particle_data>    load_particle_data();
    std::vector<index_range>      load_nucleolus_ranges();
    std::vector<nucleolus_bond>   load_nucleolus_bonds();

    // Snapshot
    void set_phase(std::string const& phase);
    void save_chromosomes(md::array_view<chromosome_range const> chroms);
    void save_positions(md::step step, md::array_view<md::point const> positions);
    void save_context(md::step step, simulation_context const& context);
    void save_contacts(md::step step, std::vector<std::array<std::uint32_t, 3>> const& contacts);

    std::vector<md::point> load_positions(md::step step);
    simulation_context     load_context(md::step step);

private:
    H5::File _store;
    std::string _phase = "unknown";
};
