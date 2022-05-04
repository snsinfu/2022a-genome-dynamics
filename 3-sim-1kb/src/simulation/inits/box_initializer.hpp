#pragma once

#include <random>
#include <vector>

#include <md.hpp>

#include "initializer.hpp"
#include "../topology.hpp"


class box_initializer : public simulation_initializer
{
public:
    struct config_type
    {
        md::scalar                    box_size    = 0;
        md::scalar                    bond_length = 0;
        std::vector<chain_assignment> chains;
    };

    explicit box_initializer(config_type const& config);

    void initialize(md::system& system, std::mt19937_64& random) const override;

private:
    config_type _config;
};
