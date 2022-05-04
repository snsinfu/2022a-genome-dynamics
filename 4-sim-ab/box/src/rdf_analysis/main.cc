#define DOCOPT_HEADER_ONLY

#include <iostream>
#include <string>

#include <docopt/docopt.h>
#include <md.hpp>

#include "analysis.hpp"


namespace
{
    char const usage[] = R"(
usage:
  rdf_analysis [options] <FILE>

  <FILE>  HDF5 trajectory file to analyze

options:
  --type <TYPE>          Particle type to select
  --steps <RANGE>        Step or range of steps (start:end) to analyze
  --bin-width <DIST>     Bin width
  --max-distance <DIST>  Max distance for analysis
  -h, --help             Print this help message and exit
)";

    analysis_config parse_options(int argc, char** argv);
}


int main(int argc, char** argv)
{
    try {
        run_analysis(parse_options(argc, argv));
    } catch (std::exception const& e) {
        std::cerr << "error: " << e.what() << '\n';
    }
}


namespace
{
    struct integer_range
    {
        long start;
        long end;
    };

    integer_range parse_range(std::string const& arg);
    md::scalar    parse_distance(std::string const& arg);


    analysis_config parse_options(int argc, char** argv)
    {
        analysis_config config;

        auto const options = docopt::docopt(usage, {argv + 1, argv + argc});

        config.filename = options.at("<FILE>").asString();

        if (auto const option = options.at("--type")) {
            config.particle_selection = option.asString();
        }

        if (auto const option = options.at("--steps")) {
            auto [ start, end ] = parse_range(option.asString());
            config.step_start = static_cast<md::step>(start);
            config.step_end = static_cast<md::step>(end);
        }

        if (auto const option = options.at("--bin-width")) {
            config.bin_width = parse_distance(option.asString());
        }

        if (auto const option = options.at("--max-distance")) {
            config.max_distance = parse_distance(option.asString());
        }

        return config;
    }


    // Parses integer range in the form `start:end`. `start` can be omitted
    // and defaults to zero.
    integer_range parse_range(std::string const& arg)
    {
        std::size_t pos;
        auto start = std::stol(arg, &pos);
        auto end = start;

        if (pos < arg.size()) {
            if (arg[pos] == ':') {
                // Range specification like 100:500
                end = std::stol(arg.substr(pos + 1));
            } else {
                throw std::invalid_argument("invalid range specification");
            }
        }

        return {start, end};
    }


    // Parses distance parameter.
    md::scalar parse_distance(std::string const& arg)
    {
        return std::stod(arg);
    }
}
