#pragma once

#include <random>

#include <md.hpp>


class simulation_initializer
{
public:
    virtual ~simulation_initializer() = default;
    virtual void initialize(md::system& system, std::mt19937_64& random) const = 0;
};
