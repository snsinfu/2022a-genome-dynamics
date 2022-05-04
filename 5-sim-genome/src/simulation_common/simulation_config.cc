#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <md.hpp>
#include <nlohmann/json.hpp>

#include "simulation_config.hpp"


namespace
{
    template<typename T>
    void load(nlohmann::json const& node, T& var);

    template<typename T>
    void store(nlohmann::json& node, T const& val);
}


void parse_simulation_config(std::string const& str, simulation_config& config)
{
    auto const json = nlohmann::json::parse(str);

    detail::foreach_simulation_config_parameter(
        config,
        [&](std::string const& name, auto& var) {
            if (auto node = json.find(name); node != json.end()) {
                load(*node, var);
            } else {
                throw std::runtime_error(name + " is not configured");
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
            // Needs c_str for suppressing a ccls error.
            store(json[name.c_str()], var);
        }
    );
    return json.dump(/*pretty=*/ true);
}


namespace
{
    // Function: load
    //
    // Loads a value of type T from a JSON node.
    //
    template<typename T>
    void load(nlohmann::json const& node, T& var)
    {
        var = node;
    }

    template<>
    void load<md::vector>(nlohmann::json const& node, md::vector& var)
    {
        std::vector<nlohmann::json> coords = node;
        var = {coords.at(0), coords.at(1), coords.at(2)};
    }

    template<>
    void load<md::point>(nlohmann::json const& node, md::point& var)
    {
        std::vector<nlohmann::json> coords = node;
        var = {coords.at(0), coords.at(1), coords.at(2)};
    }


    // Function: store
    //
    // Stores a value of type T to a JSON node.
    //
    template<typename T>
    void store(nlohmann::json& node, T const& val)
    {
        node = val;
    }

    template<>
    void store<md::vector>(nlohmann::json& node, md::vector const& val)
    {
        node = std::vector<nlohmann::json>{val.x, val.y, val.z};
    }

    template<>
    void store<md::point>(nlohmann::json& node, md::point const& val)
    {
        node = std::vector<nlohmann::json>{val.x, val.y, val.z};
    }
}
