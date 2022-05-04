#pragma once

// This module provides simulation_config struct and its (de)serialization
// functions.

#include <cstdint>
#include <string>

#include <md.hpp>


// Structure: simulation_config
//
// Simple aggregate of simulation parameters.
//
struct simulation_config
{
#define X(var, T, value)    T var = value;
#define V(...)              {__VA_ARGS__}
#include "../config_entries.inc"
#undef X
#undef V
};


namespace detail
{
    // Function: foreach_simulation_config_parameter
    //
    // Calls op(name, var) for each parameter variable in config. Config can be
    // simulation_config or const-qualified one. This function is used to
    // implement `load/dump_simulation_config` functions.
    //
    template<typename Config, typename Op>
    void foreach_simulation_config_parameter(Config& config, Op op)
    {
#define X(var, T, value)    op(#var, config.var);
#define V(...)
#include "../config_entries.inc"
#undef X
#undef V
    }
}


// Function: parse_simulation_config
//
// Loads parameter values from JSON string. Config entries are untouched if not
// listed in the JSON.
//
void parse_simulation_config(std::string const& str, simulation_config& config);


// Function: dump_simulation_config
//
// Dumps parameter values as a JSON string.
//
std::string dump_simulation_config(const simulation_config& config);
