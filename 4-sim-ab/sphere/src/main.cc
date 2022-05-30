#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#define DOCOPT_HEADER_ONLY
#include <docopt/docopt.h>

#include "simulation_config.hpp"
#include "simulation_driver.hpp"


namespace
{
    using options_map = std::map<std::string, docopt::value>;

    simulation_config make_config(options_map const& options);
}

char const program_usage[] = R"(
usage:
  simulation [-s <seed>] <config> <out>

  <config>    input JSON configuration file
  <out>       output HDF5 trajectory file

options:
  -s <seed>   specify random seed
  -h, --help  print this usage message
)";


int main(int argc, char** argv)
{
    auto const options = docopt::docopt(program_usage, {argv + 1, argv + argc});
    auto const config = make_config(options);

    simulation_driver sim{config};
    sim.run();
}


namespace
{
    // Creates `simulation_config` from program options.
    simulation_config make_config(options_map const& options)
    {
        simulation_config config;

        if (std::ifstream config_file{options.at("<config>").asString()}) {
            load_simulation_config(config_file, config);
        } else {
            throw std::runtime_error("cannot open config file");
        }

        if (auto seed_option = options.at("-s")) {
            config.seed = static_cast<std::uint64_t>(seed_option.asLong());
        }
        config.output = options.at("<out>").asString();

        return config;
    }
}
