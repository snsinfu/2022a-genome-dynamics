#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include <md.hpp>


struct sampling_config
{
    md::scalar    temperature          = 1;
    md::scalar    timestep             = 0;
    md::step      steps                = 0;
    md::step      clear_loops_at       = 0;
    bool          loop_preloading      = false;
    md::step      loop_update_interval = 1;
    md::step      glue_update_interval = 1;
    md::step      logging_interval     = 1;
    md::step      sampling_interval    = 1;
    std::uint64_t random_seed          = 0;
    std::string   output_filename;
};


struct chain_type_config
{
    md::scalar box_size            = 1;
    md::scalar initial_bond_length = 0;
    md::scalar repulsive_diameter  = 0;
    md::scalar repulsive_energy    = 0;
    md::scalar attractive_diameter = 0;
    md::scalar attractive_energy   = 0;
    md::scalar bond_length         = 0;
    md::scalar bond_spring         = 0;
    md::scalar bending_energy      = 0;
    md::scalar monomer_mobility    = 1;
};


struct loop_type_config
{
    md::scalar                bond_spring              = 0;
    md::scalar                forward_speed            = 0;
    md::scalar                backward_speed           = 0;
    md::scalar                loading_rate_density     = 0;
    md::scalar                unloading_rate           = 0;
    md::scalar                convergent_detachability = 1;
    md::scalar                roadblock_attachability  = 1;
    std::optional<md::scalar> crossing_rate;
    std::optional<md::index>  max_loops;
};


struct glue_type_config
{
    md::index  max_glues           = 0;
    md::scalar glue_energy         = 0;
    md::scalar glue_distance       = 0;
    md::scalar glue_binding_rate   = 0;
    md::scalar glue_unbinding_rate = 0;
};


struct block_config
{
    md::index                 start = 0;
    md::index                 end   = 0;
    std::optional<md::scalar> bending_energy;
};


struct chain_config
{
    md::index                 length = 0;
    std::vector<md::index>    forward_boundaries;
    std::vector<md::index>    backward_boundaries;
    std::vector<md::index>    roadblocks;
    std::vector<md::index>    loaded_loops;
    std::vector<block_config> blocks;
};


struct simulation_config
{
    sampling_config           sampling;
    chain_type_config         chain;
    loop_type_config          loop;
    glue_type_config          glue;
    std::vector<chain_config> chains;
    std::string               config_text;
};


/** Parses JSON representation of `simulation_config` structure. */
simulation_config parse_simulation_config(std::string const& text);


/** Parses JSON representation of an array of `chain_config` structure. */
std::vector<chain_config> parse_chains_config(std::string const& text);


/** Formats `simulation_config` structure as a JSON string. */
std::string format_simulation_config(simulation_config const& config);
