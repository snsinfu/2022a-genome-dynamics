#include <md.hpp>

#include "loop_forcefield.hpp"
#include "../loops/loop_simulator.hpp"


loop_forcefield::loop_forcefield(
    std::shared_ptr<loop_simulator> const& simulator,
    md::scalar                             stiffness,
    md::scalar                             length
)
    : _simulator{simulator}
    , _potential{stiffness, length}
{
}


md::scalar
loop_forcefield::compute_energy(md::system const& system)
{
    auto const positions = system.view_positions();
    md::scalar energy = 0;

    for (auto const& loop : *_simulator) {
        if (!loop.id) {
            continue;
        }
        auto const r_ij = positions[loop.start] - positions[loop.end];
        energy += _potential.evaluate_energy(r_ij);
    }

    return energy;
}


void
loop_forcefield::compute_force(md::system const& system, md::array_view<md::vector> forces)
{
    auto const positions = system.view_positions();

    for (auto const& loop : *_simulator) {
        if (!loop.id) {
            continue;
        }
        auto const r_ij = positions[loop.start] - positions[loop.end];
        auto const f_ij = _potential.evaluate_force(r_ij);
        forces[loop.start] += f_ij;
        forces[loop.end] -= f_ij;
    }
}
