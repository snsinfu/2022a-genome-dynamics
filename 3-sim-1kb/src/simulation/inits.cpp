#include <stdexcept>
#include <memory>

#include "config.hpp"
#include "inits.hpp"
#include "inits/box_initializer.hpp"
#include "topology.hpp"


std::unique_ptr<simulation_initializer>
make_initializer(simulation_config const& config)
{
    auto const chains = make_chain_assignments(config);

    return std::make_unique<box_initializer>(
        box_initializer::config_type{
            .box_size    = config.chain.box_size,
            .bond_length = config.chain.initial_bond_length,
            .chains      = chains,
        }
    );
}
