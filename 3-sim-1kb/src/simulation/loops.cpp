#include <cstddef>
#include <memory>
#include <vector>

#include "config.hpp"
#include "loops.hpp"
#include "loops/basic_loop_simulator.hpp"
#include "topology.hpp"


std::shared_ptr<loop_simulator>
make_loop_simulator(simulation_config const& config)
{
    auto const chains = make_chain_assignments(config);

    // We simulate 1D loop dynamics on a concatenated virtual chain. Loop
    // factors do not hop across boundaries. So, as long as boundary elements
    // are correctly assigned on the virtual chain, the 1D loop dynamics
    // would be correct.
    std::size_t virtual_length = 0;
    for (auto const& chain : chains) {
        virtual_length += chain.config.length;
    }

    std::size_t max_loops = 0;
    if (config.loop.max_loops) {
        max_loops = *config.loop.max_loops;
    } else {
        // Deduce from the number of initially loaded loops.
        for (auto const& chain : chains) {
            max_loops += chain.config.loaded_loops.size();
        }
    }

    basic_loop_simulator loops{{
        .chain_length = virtual_length,
        .max_loops    = max_loops,
    }};

    loops.set_forward_speed(config.loop.forward_speed);
    loops.set_backward_speed(config.loop.backward_speed);
    loops.set_loading_rate(config.loop.loading_rate_density * double(virtual_length));
    loops.set_unloading_rate(config.loop.unloading_rate);

    if (config.loop.crossing_rate) {
        loops.set_crossing_rate(*config.loop.crossing_rate);
    }

    for (auto const& chain : chains) {
        // Add forward boundaries. Loop factors slow down at the right neighbor
        // of each forward boundary.
        for (auto const& pos : chain.config.forward_boundaries) {
            auto const virtual_pos = chain.start + pos;
            loops.add_boundary(virtual_pos);
            if (virtual_pos + 1 < chain.end) {
                loops.set_site_detachability(
                    virtual_pos + 1,
                    config.loop.convergent_detachability
                );
            }
        }

        // Add backward boundaries. Loop factors slow down at the left neighbor
        // of each backward boundary.
        for (auto const& pos : chain.config.backward_boundaries) {
            auto const virtual_pos = chain.start + pos;
            loops.add_boundary(virtual_pos);
            if (virtual_pos >= chain.start + 1) {
                loops.set_site_detachability(
                    virtual_pos - 1,
                    config.loop.convergent_detachability
                );
            }
        }

        // Add roadblocks. Loop factors enter these sites in lower probability.
        for (auto const& pos : chain.config.roadblocks) {
            auto const virtual_pos = chain.start + pos;
            loops.set_site_attachability(virtual_pos, config.loop.roadblock_attachability);
        }
    }

    // Load initial, "handcuff" loops.
    for (auto const& chain : chains) {
        for (auto const& pos : chain.config.loaded_loops) {
            auto const virtual_pos = chain.start + pos;
            loops.load_loop(virtual_pos);
        }
    }

    return std::make_shared<basic_loop_simulator>(loops);
}
