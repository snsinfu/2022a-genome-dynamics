#include <vector>

#include "topology.hpp"


std::vector<chain_assignment> make_chain_assignments(simulation_config const& config)
{
    std::vector<chain_assignment> assignments;
    md::index offset = 0;

    for (auto const& chain : config.chains) {
        assignments.push_back({
            .start  = offset,
            .end    = offset + chain.length,
            .config = chain,
        });
        offset += chain.length;
    }

    return assignments;
}
