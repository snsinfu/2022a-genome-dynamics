#define DOCOPT_HEADER_ONLY

#include <exception>
#include <iostream>
#include <string>

#include <docopt/docopt.h>
#include <md.hpp>

#include "analysis.hpp"


static char const usage[] = R"(
usage:
  rdf_analysis_hetero [options] <FILE>

  <FILE>  HDF5 trajectory file to analyze

options:
  --type <TYPE>          Type of center particles used in RDF analysis [default: A]
  --bin-width <DIST>     Bin width
  --max-distance <DIST>  Max distance for analysis
  -h, --help             Print this help message and exit
)";

static analysis_config parse_options(int argc, char** argv);


int main(int argc, char** argv)
{
    try {
        run_analysis(parse_options(argc, argv));
    } catch (std::exception const& e) {
        std::cerr << "error: " << e.what() << '\n';
    }
}


static analysis_config parse_options(int argc, char** argv)
{
    analysis_config config;

    auto const options = docopt::docopt(usage, {argv + 1, argv + argc});

    config.filename = options.at("<FILE>").asString();

    if (auto const option = options.at("--type")) {
        config.center_type = option.asString();
    }

    if (auto const option = options.at("--bin-width")) {
        config.bin_width = std::stod(option.asString());
    }

    if (auto const option = options.at("--max-distance")) {
        config.max_distance = std::stod(option.asString());
    }

    return config;
}
