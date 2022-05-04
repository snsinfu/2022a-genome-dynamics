#include <random>
#include <vector>

#include <md.hpp>

#include "../topology.hpp"


inline void
generate_random_walks(
    md::array_view<md::point>            positions,
    std::vector<chain_assignment> const& chains,
    std::vector<md::point> const&        centroids,
    md::scalar                           bond_length,
    std::mt19937_64&                     random
)
{
    auto direction = [](auto& random) {
        std::normal_distribution<md::scalar> normal;
        return md::normalize(md::vector{normal(random), normal(random), normal(random)});
    };

    // Generate a random walk path.
    md::point walk;

    for (auto& pos : positions) {
        pos = walk;
        walk += bond_length * direction(random);
    }

    // Correct the centroid of each chain to the desired position.
    for (std::size_t chain_index = 0; chain_index < chains.size(); chain_index++) {
        auto const& chain = chains[chain_index];
        auto const centroid = centroids[chain_index];
        auto const chain_length = chain.end - chain.start;

        md::vector offset;
        for (auto const& pos : positions.subview(chain.start, chain_length)) {
            offset += pos - centroid;
        }
        offset /= md::scalar(chain_length);

        for (auto& pos : positions.subview(chain.start, chain_length)) {
            pos -= offset;
        }
    }
}
