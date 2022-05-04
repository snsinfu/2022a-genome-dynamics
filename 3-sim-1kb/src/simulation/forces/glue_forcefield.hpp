#pragma once

#include <memory>

#include <md.hpp>

#include "../glues/glue_simulator.hpp"


class glue_forcefield : public md::forcefield
{
public:
    explicit glue_forcefield(
        std::shared_ptr<glue_simulator> const& simulator,
        md::periodic_box const& box,
        md::scalar energy,
        md::scalar cutoff_distance
    );

    md::scalar compute_energy(md::system const& system) override;
    void       compute_force(md::system const& system, md::array_view<md::vector> forces) override;

private:
    std::shared_ptr<glue_simulator> _simulator;
    md::periodic_box                _box;
    md::softcore_potential<8, 3>    _potential;
};