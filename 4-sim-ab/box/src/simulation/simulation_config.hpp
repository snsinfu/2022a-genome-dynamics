#pragma once

#include <cstdint>
#include <istream>
#include <string>

#include <md.hpp>


// Enumerate (de)serielizable parameters here. X macro idiom.
#define X_CONFIG_JSON_PARAMETERS              \
    X(md::scalar,    a_core_diameter,   0.30) \
    X(md::scalar,    b_core_diameter,   0.24) \
    X(md::scalar,    a_core_repulsion,  2.0 ) \
    X(md::scalar,    b_core_repulsion,  2.0 ) \
    X(md::scalar,    bond_spring,       70  ) \
    X(md::scalar,    mobility,          1.0 ) \
    X(md::scalar,    box_size,          1.0 ) \
    X(md::scalar,    init_bond_length,  0.0 ) \
    X(md::scalar,    temperature,       1.0 ) \
    X(md::scalar,    timestep,          1e-5) \
    X(md::step,      steps,             1000) \
    X(md::step,      logging_interval,  1000) \
    X(md::step,      sampling_interval, 1000) \
    X(std::string,   beads_filename,    ""  ) \
    X(std::uint64_t, seed,              0   )


// Simple aggregate of simulation parameters.
struct simulation_config
{
    std::string output;

#define X(T, var, init) T var = init;
    X_CONFIG_JSON_PARAMETERS
#undef X
};


namespace detail
{
    // Calls op(name, var) for each parameter variable in config. Config can be
    // simulation_config or const-qualified one.
    template<typename Config, typename Op>
    void foreach_simulation_config_parameter(Config& config, Op op)
    {
#define X(T, var, init) op(#var, config.var);
        X_CONFIG_JSON_PARAMETERS
#undef X
    }
}


// Loads parameter values from JSON input. Config entries are untouched if they
// are not listed in the JSON.
void load_simulation_config(std::istream& in, simulation_config& config);


// Dumps parameter values as a JSON string. The `output` member is not dumped.
std::string dump_simulation_config(simulation_config const& config);


#undef X_CONFIG_JSON_PARAMETERS
