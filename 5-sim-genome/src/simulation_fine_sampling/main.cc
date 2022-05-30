#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <md.hpp>

#include "../simulation_common/simulation_store.hpp"

#include "simulation_driver.hpp"


int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: simulation_fine_sampling <trajectory>\n";
        return 1;
    }

    simulation_store store{argv[1]};
    simulation_driver driver{store};
    driver.run();

    return 0;
}
