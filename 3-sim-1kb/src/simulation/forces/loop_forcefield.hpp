#pragma once

#include <memory>

#include <md.hpp>

#include "../loops/loop_simulator.hpp"


class loop_forcefield : public md::forcefield
{
public:
    explicit loop_forcefield(
        std::shared_ptr<loop_simulator> const& simulator,
        md::scalar stiffness,
        md::scalar length
    );

    md::scalar compute_energy(md::system const& system) override;
    void compute_force(md::system const& system, md::array_view<md::vector> forces) override;

private:
    std::shared_ptr<loop_simulator> _simulator;
    md::spring_potential            _potential;
};
