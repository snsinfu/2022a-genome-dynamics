#include <string>
#include <vector>

#include <jsoncons/json.hpp>

#include "config.hpp"


JSONCONS_N_MEMBER_TRAITS(
    sampling_config,

    // Required fields
    3,
    temperature,
    timestep,
    steps,

    // Optional fields
    clear_loops_at,
    loop_preloading,
    loop_update_interval,
    glue_update_interval,
    logging_interval,
    sampling_interval,
    random_seed,
    output_filename
)


JSONCONS_N_MEMBER_TRAITS(
    chain_type_config,

    // Required fields
    0,

    // Optional fields
    box_size,
    initial_bond_length,
    repulsive_diameter,
    repulsive_energy,
    attractive_diameter,
    attractive_energy,
    bond_length,
    bond_spring,
    bending_energy,
    monomer_mobility
)


JSONCONS_N_MEMBER_TRAITS(
    loop_type_config,

    // Required fields
    3,
    bond_spring,
    forward_speed,
    backward_speed,

    // Optional fields
    loading_rate_density,
    unloading_rate,
    convergent_detachability,
    roadblock_attachability,
    crossing_rate,
    max_loops
)


JSONCONS_ALL_MEMBER_TRAITS(
    glue_type_config,

    // Required fields
    max_glues,
    glue_energy,
    glue_distance,
    glue_binding_rate,
    glue_unbinding_rate
)


JSONCONS_N_MEMBER_TRAITS(
    chain_config,

    // Required fields
    1,
    length,

    // Optional fields
    forward_boundaries,
    backward_boundaries,
    roadblocks,
    loaded_loops,
    blocks
)


JSONCONS_N_MEMBER_TRAITS(
    block_config,

    // Required fields
    2,
    start,
    end,

    // Optional fields
    bending_energy
)


JSONCONS_N_MEMBER_TRAITS(
    simulation_config,

    // Required fields
    1,
    sampling,

    // Optional fields
    loop,
    chain,
    glue,
    chains
)


simulation_config parse_simulation_config(std::string const& text)
{
    auto config = jsoncons::decode_json<simulation_config>(text);
    config.config_text = text;
    return config;
}


std::vector<chain_config> parse_chains_config(std::string const& text)
{
    return jsoncons::decode_json<std::vector<chain_config>>(text);
}


std::string format_simulation_config(simulation_config const& config)
{
    std::string text;
    jsoncons::encode_json(config, text);
    return text;
}
