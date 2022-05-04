#include <algorithm>
#include <cmath>
#include <iterator>

#include "function_output_iterator.hpp"
#include "glue_simulator.hpp"
#include "reservoir_sampler.hpp"


glue_simulator::glue_simulator(config_type const& config)
    : _config{config}
    , _searcher{config.box, config.max_distance}
{
}

std::size_t
glue_simulator::size() const
{
    return _glued_pairs.size();
}

glue_simulator::iterator
glue_simulator::begin() const
{
    return _glued_pairs.begin();
}

glue_simulator::iterator
glue_simulator::end() const
{
    return _glued_pairs.end();
}

void
glue_simulator::update(double timestep, md::array_view<md::point const> positions, random_engine& random)
{
    if (_config.max_glues == 0) {
        return;
    }

    _searcher.set_points(positions);
    unbond(timestep, positions, random);
    bond(timestep, positions, random);
}

void
glue_simulator::unbond(double timestep, md::array_view<md::point const> positions, random_engine& random)
{
    for (auto next = _glued_pairs.cbegin(); next != _glued_pairs.cend(); ) {
        // Need to fetch next node before (potentially) removing the current node.
        auto const pair = next++;

        auto const r_ij = _config.box.shortest_displacement(positions[pair->i], positions[pair->j]);
        std::bernoulli_distribution unbinding{-std::expm1(-_config.unbinding_rate * timestep)};

        if (r_ij.norm() > _config.max_distance || unbinding(random)) {
            _glued_pairs.erase(pair);
        }
    }
}

void
glue_simulator::bond(double timestep, md::array_view<md::point const>, random_engine& random)
{
    reservoir_sampler<glue_pair> reservoir{_config.max_glues - _glued_pairs.size()};

    _searcher.search(function_output_iterator([&](auto const& ij) {
        glue_pair const pair = {
            .i = std::uint32_t(ij.first),
            .j = std::uint32_t(ij.second),
        };
        std::bernoulli_distribution binding{-std::expm1(-_config.binding_rate * timestep)};

        if (!_glued_pairs.contains(pair) && binding(random)) {
            reservoir.feed(pair, random);
        }
    }));

    std::copy(reservoir.begin(), reservoir.end(), std::inserter(_glued_pairs, _glued_pairs.end()));
}
