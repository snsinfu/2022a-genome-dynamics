#include <md.hpp>

#include "../simulation_common/particle_data.hpp"

#include "simulation_driver.hpp"


void simulation_driver::setup_particles()
{
    auto const particles = _store.load_particle_data();
    auto const chromosomes = _store.load_chromosomes();
    auto const nucleolus_ranges = _store.load_nucleolus_ranges();

    // Just copy particle data from the simulation metadata.
    _system.add_attribute(particle_data_attribute);

    for (md::index i = 0; i < particles.size(); i++) {
        auto part = _system.add_particle();
        part.view(particle_data_attribute) = particles[i];
    }

    // Mobility varies between chromatin and nucleolar particle.
    auto mobilities = _system.view_mobilities();

    for (auto const& chrom : chromosomes) {
        for (md::index i = chrom.start; i < chrom.end; i++) {
            mobilities[i] = _config.chromatin_mobility;
        }
    }

    for (auto const& nucleo : nucleolus_ranges) {
        for (md::index i = nucleo.begin; i < nucleo.end; i++) {
            mobilities[i] = _config.nucleolus_mobility;
        }
    }
}
