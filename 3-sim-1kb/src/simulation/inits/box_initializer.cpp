#include <random>
#include <vector>

#include <md.hpp>

#include "box_initializer.hpp"
#include "utils.hpp"
#include "../topology.hpp"


box_initializer::box_initializer(config_type const& config)
    : _config{config}
{
}


void
box_initializer::initialize(md::system& system, std::mt19937_64& random) const
{
    // Place chains at random positions in the box, and initialize each chain
    // as random walk paths.

    std::vector<md::point> centroids;

    for (auto const& chain : _config.chains) {
        std::uniform_real_distribution<md::scalar> coord{0, _config.box_size};
        centroids.push_back({coord(random), coord(random), coord(random)});
        (void) chain;
    }

    generate_random_walks(
        system.view_positions(),
        _config.chains,
        centroids,
        _config.bond_length,
        random
    );
}
