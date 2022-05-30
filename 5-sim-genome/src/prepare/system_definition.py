from dataclasses import dataclass
from typing import List, Tuple


# Particle type codes.
TYPE_A = 1
TYPE_B = 2
TYPE_U = 3
TYPE_CEN = 4
TYPE_ANOR = 5
TYPE_BNOR = 6
TYPE_NUC = 7

# Mapping among type code, genome tag and descriptive name. The order matters.
# Preceding type is taken for multi-tagged region.
TYPE_CODEBOOK = [
    (TYPE_ANOR, "anor", "active_NOR"),
    (TYPE_BNOR, "bnor", "silent_NOR"),
    (TYPE_CEN, "cen", "centromere"),
    (TYPE_A, "A", "A"),
    (TYPE_B, "B", "B"),
    (TYPE_U, "u", "u"),
    (TYPE_NUC, None, "nucleolus"),
]


@dataclass
class Particle:
    type: int
    A: float
    B: float
    index: int


@dataclass
class Span:
    name: str
    start: int
    end: int


@dataclass
class ChromChain(Span):
    cen_start: int
    cen_end: int


@dataclass
class SystemDefinition:
    particles: List[Particle]
    chromatin_chains: List[ChromChain]
    nucleolus_spans: List[Span]
    nucleolus_bonds: List[Tuple[int, int]]


def make_system_definition(genome, config):
    """
    Create a new SystemDefinition object from given genome table.
    """
    system = SystemDefinition([], [], [], [])
    _add_chromatin_particles(system, genome, config)
    _add_nucleolus_particles(system, genome, config)
    return system


def _add_chromatin_particles(system, genome, config):
    """
    Define and add chromatin particles to the system.
    """
    particles = system.particles
    chromatin_chains = system.chromatin_chains

    for name, chain in genome.groupby(genome.chain, sort=False):
        start = len(particles)
        cens = []

        for _, bead in chain.iterrows():
            type_code = _compute_type_code(bead)
            index = len(particles)
            particles.append(Particle(type_code, bead.A, bead.B, index))

            if type_code == TYPE_CEN:
                cens.append(index)

        end = len(particles)

        if cens:
            cen_start = cens[0]
            cen_end = cens[-1]
        else:
            cen_start = cen_end = 0
        chromatin_chains.append(ChromChain(name, start, end, cen_start, cen_end))


def _add_nucleolus_particles(system, genome, config):
    """
    Define and add nucleolar particles to the system.
    """
    particles = system.particles
    chromatin_chains = system.chromatin_chains
    nucleolus_spans = system.nucleolus_spans
    nucleolus_bonds = system.nucleolus_bonds

    for chain in chromatin_chains:
        start = len(particles)

        for nor_index in range(chain.start, chain.end):
            # Nucleolar particles are only attached to active NORs.
            if particles[nor_index].type != TYPE_ANOR:
                continue

            for _ in range(config["nucleolus_sidebeads"]):
                nuc_index = len(particles)
                particles.append(
                    Particle(
                        TYPE_NUC,
                        config["nucleolus_a_factor"],
                        config["nucleolus_b_factor"],
                        nuc_index,
                    )
                )
                nucleolus_bonds.append((nor_index, nuc_index))

        end = len(particles)

        if start != end:
            nucleolus_spans.append(Span(chain.name, start, end))


def _compute_type_code(bead):
    """
    Determine particle type from chromatin data (a row in a genome table).
    """
    tags = bead.tags.split(",")

    for code, tag, _ in TYPE_CODEBOOK:
        if tag in tags:
            return code

    assert False
