#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.hpp>

#include "config.hpp"
#include "simulation.hpp"


/** Program mode specified by user. */
enum class program_mode
{
    simulation, // Run simulation
    help,       // Show command-line usage and exit
};


/** Command-line arguments and option flags. */
struct program_options
{
    program_mode                 mode;
    std::string                  config_filename;
    std::optional<std::string>   output_filename;
    std::optional<std::string>   chains_filename;
    std::optional<std::uint64_t> random_seed;
};


static void                      show_usage();
static program_options           parse_options(int argc, char** argv);
static simulation_config         make_config(program_options const& options);
static simulation_config         load_config(std::string const& filename);
static std::vector<chain_config> load_chains(std::string const& filename);
static std::string               load_text(std::string const& filename);


int main(int argc, char** argv)
{
    try {
        auto const options = parse_options(argc, argv);

        switch (options.mode) {
        case program_mode::simulation:
            simulation{make_config(options)}.run();
            break;

        case program_mode::help:
            show_usage();
            break;
        }

        return 0;
    } catch (std::exception const& err) {
        std::cerr << "error: " << err.what() << '\n';
        return 1;
    }
}


void show_usage()
{
    std::cerr <<
        "Loop formation simulator\n"
        "usage: main [-hCos] <config>\n"
        "\n"
        "  <config>     JSON file specifying simulation parameters\n"
        "\n"
        "options:\n"
        "  -C <config>  override chain definitions (config 'chains' key) by additional JSON file\n"
        "  -o <output>  override output HDF5 filename (config 'output_filename' key)\n"
        "  -s <seed>    override random seed (config 'random_seed' key)\n"
        "  -h           print this usage message and exit\n"
        "\n";
}


/** Parses command-line arguments and build a `program_options` structure. */
program_options parse_options(int argc, char** argv)
{
    program_options options;
    options.mode = program_mode::simulation;

    // Optional arguments

    cxx::getopt getopt;

    for (int opt; (opt = getopt(argc, argv, "C:o:s:h")) != -1; ) {
        switch (opt) {
        case 'C':
            options.chains_filename = getopt.optarg;
            break;

        case 'o':
            options.output_filename = getopt.optarg;
            break;

        case 's':
            options.random_seed = std::stoull(getopt.optarg);
            break;

        case 'h':
            options.mode = program_mode::help;
            // Do not parse other arguments if help is requested.
            return options;

        case '?':
            throw std::runtime_error{"bad option"};
        }
    }

    argc -= getopt.optind;
    argv += getopt.optind;

    // Positional arguments

    if (argc != 1) {
        throw std::runtime_error{"config file is not specified"};
    }

    options.config_filename = argv[0];

    return options;
}


/** Builds simulation configuration from given program options. */
simulation_config make_config(program_options const& options)
{
    auto config = load_config(options.config_filename);

    if (options.chains_filename) {
        config.chains = load_chains(*options.chains_filename);
    }

    if (options.output_filename) {
        config.sampling.output_filename = *options.output_filename;
    }

    if (options.random_seed) {
        config.sampling.random_seed = *options.random_seed;
    }

    return config;
}


/** Loads simulation configuration from a JSON file. */
simulation_config load_config(std::string const& filename)
{
    auto const text = load_text(filename);
    try {
        return parse_simulation_config(text);
    } catch (std::exception const& err) {
        throw std::runtime_error{"failed to parse config file - " + std::string{err.what()}};
    }
}


/** Loads simulation configuration from a JSON file. */
std::vector<chain_config> load_chains(std::string const& filename)
{
    auto const text = load_text(filename);
    try {
        return parse_chains_config(text);
    } catch (std::exception const& err) {
        throw std::runtime_error{"failed to parse chains config file - " + std::string{err.what()}};
    }
}


/** Loads text string from a file. */
std::string load_text(std::string const& filename)
{
    std::ifstream file{filename};
    std::string text;
    if (!std::getline(file, text, '\0')) {
        throw std::runtime_error{"failed to load config file"};
    }
    return text;
}
