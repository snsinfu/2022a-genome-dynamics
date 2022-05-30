#include "glue_forcefield.hpp"


glue_forcefield::glue_forcefield(
    std::shared_ptr<glue_simulator> const& simulator,
    md::periodic_box const& box,
    md::scalar energy,
    md::scalar cutoff_distance
)
    : _simulator{simulator}
    , _box{box}
    , _potential{.energy = -energy, .diameter = cutoff_distance}
{
}

md::scalar
glue_forcefield::compute_energy(md::system const& system)
{
    auto const positions = system.view_positions();

    md::scalar total_energy = 0;

    for (auto const& pair : *_simulator) {
        auto const r_ij = _box.shortest_displacement(positions[pair.i], positions[pair.j]);
        auto const e_ij = _potential.evaluate_energy(r_ij);
        total_energy += e_ij;
    }

    return total_energy;
}

void
glue_forcefield::compute_force(md::system const& system, md::array_view<md::vector> forces)
{
    auto const positions = system.view_positions();

    for (auto const& pair : *_simulator) {
        auto const r_ij = _box.shortest_displacement(positions[pair.i], positions[pair.j]);
        auto const f_ij = _potential.evaluate_force(r_ij);

        forces[pair.i] += f_ij;
        forces[pair.j] -= f_ij;
    }
}