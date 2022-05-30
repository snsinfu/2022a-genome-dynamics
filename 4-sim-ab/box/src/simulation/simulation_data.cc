#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "simulation_data.hpp"


std::vector<bead_data> load_beads_data(std::string const& filename)
{
    std::ifstream file{filename};

    std::string const header_expect = "chain\tA\tB";
    std::string header;
    std::getline(file, header);

    if (header != header_expect) {
        throw std::runtime_error("unexpected beads data header");
    }

    std::vector<bead_data> beads;

    for (std::string line; std::getline(file, line); ) {
        std::istringstream record{line};

        std::string chain;
        double a_factor;
        double b_factor;

        record
            >> chain
            >> a_factor
            >> b_factor;

        beads.push_back({
            .chain    = chain,
            .a_factor = a_factor,
            .b_factor = b_factor,
        });
    }

    return beads;
}
