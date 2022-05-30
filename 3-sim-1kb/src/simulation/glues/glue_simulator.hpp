#pragma once

#include <random>
#include <unordered_set>
#include <vector>

#include <md.hpp>


struct glue_pair
{
    std::uint32_t i = 0;
    std::uint32_t j = 0;
};

inline bool operator==(glue_pair const& a, glue_pair const& b)
{
    return a.i == b.i && a.j == b.j;
}

inline bool operator!=(glue_pair const& a, glue_pair const& b)
{
    return !(a == b);
}

namespace std
{
    template<>
    struct hash<glue_pair>
    {
        std::size_t operator()(glue_pair const& pair) const noexcept
        {
            std::size_t const a = pair.i & pair.j;
            std::size_t const b = pair.i ^ pair.j;
            return (a << 32) | b;
        }
    };
}

class glue_simulator
{
public:
    using random_engine = std::mt19937_64;
    using iterator = std::unordered_set<glue_pair>::const_iterator;

    struct config_type
    {
        std::size_t      max_glues      = 0;
        double           max_distance   = 0;
        double           binding_rate   = 0;
        double           unbinding_rate = 0;
        md::periodic_box box;
    };

    explicit glue_simulator(config_type const& config);

    std::size_t size() const;
    iterator begin() const;
    iterator end() const;

    void update(
        double                          timestep,
        md::array_view<md::point const> positions,
        random_engine&                  random
    );

private:
    void unbond(
        double                          timestep,
        md::array_view<md::point const> positions,
        random_engine&                  random
    );

    void bond(
        double                          timestep,
        md::array_view<md::point const> positions,
        random_engine&                  random
    );

private:
    config_type                             _config;
    md::neighbor_searcher<md::periodic_box> _searcher;
    std::vector<md::point>                  _positions;
    std::unordered_set<glue_pair>           _glued_pairs;
};
