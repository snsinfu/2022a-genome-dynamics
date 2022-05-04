#include <iostream>
#include <istream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "simulation_config.hpp"


void load_simulation_config(std::istream& in, simulation_config& config)
{
    auto const& json = nlohmann::json::parse(in);

    detail::foreach_simulation_config_parameter(
        config,
        [&](std::string const& name, auto& var) {
            if (auto node = json.find(name); node != json.end()) {
                var = *node;
            }
        }
    );
}


std::string dump_simulation_config(simulation_config const& config)
{
    nlohmann::json json;

    detail::foreach_simulation_config_parameter(
        config,
        [&](std::string const& name, auto const& var) {
            json[name] = var;
        }
    );
    return json.dump(/*pretty=*/ true);
}
