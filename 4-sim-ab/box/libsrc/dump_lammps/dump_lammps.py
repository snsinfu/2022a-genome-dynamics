import json
import os

import h5py
import numpy as np
import pandas as pd


PARTICLE_TYPE_A = 1
PARTICLE_TYPE_B = 2
PARTICLE_TYPE_U = 3


def run(trajfile, datafile, dumpfile, *, step_range=None):
    with h5py.File(trajfile, "r") as store:
        config = json.loads(store["metadata/config"][()])
        snapshots = store["snapshots"]

        beads_path = os.path.join(os.path.dirname(trajfile), config["beads_filename"],)
        beads = pd.read_csv(beads_path, sep="\t")

        with open(datafile, "w") as data:
            create_data(snapshots, data, config, beads)

        with open(dumpfile, "w") as dump:
            create_dump(snapshots, dump, config, beads, step_range)


def create_data(snapshots, data, config, beads):
    step = snapshots[".steps"][0]
    points = snapshots[step]["positions"][:]

    # Definition of the chains
    chain_ranges = determine_chain_ranges(beads)
    particle_types = determine_particle_types(beads)

    n_atoms = points.shape[0]
    n_bonds = sum(end - beg - 1 for beg, end in chain_ranges)
    n_atom_types = particle_types.max() - particle_types.min() + 1
    n_bond_types = 1

    # Definition of the unit cell
    x_lo = y_lo = z_lo = 0
    x_hi = y_hi = z_hi = config["box_size"]

    data.write("Simulation\n")
    data.write("\n")
    data.write(f"{n_atoms} atoms\n")
    data.write(f"{n_bonds} bonds\n")
    data.write("\n")
    data.write(f"{n_atom_types} atom types\n")
    data.write(f"{n_bond_types} bond types\n")
    data.write("\n")
    data.write(f"{x_lo} {x_hi} xlo xhi\n")
    data.write(f"{y_lo} {y_hi} ylo yhi\n")
    data.write(f"{z_lo} {z_hi} zlo zhi\n")
    data.write("\n")
    data.write("Atoms\n")
    data.write("\n")

    for i, point in enumerate(points):
        atom_id = i + 1  # Must start from 1
        atom_type = particle_types[i]
        x, y, z = point
        data.write(f"{atom_id} {atom_type} {x:g} {y:g} {z:g}\n")

    data.write("\n")
    data.write("Bonds\n")
    data.write("\n")

    bond_id = 1  # Must start from 1
    bond_type = 1
    for beg, end in chain_ranges:
        for i in range(beg, end - 1):
            j = i + 1  # Bond with the next particle
            atom_i = i + 1
            atom_j = j + 1
            data.write(f"{bond_id} {bond_type} {atom_i} {atom_j}\n")
            bond_id += 1


def create_dump(snapshots, dump, config, beads, step_range=None):
    # Definition of the chains
    chain_ranges = determine_chain_ranges(beads)
    particle_types = determine_particle_types(beads)

    # Definition of the unit cell
    x_lo = y_lo = z_lo = 0
    x_hi = y_hi = z_hi = config["box_size"]

    for step in snapshots[".steps"]:
        if step_range is not None:
            if int(step) < step_range[0]:
                continue
            if int(step) > step_range[1]:
                continue

        points = snapshots[step]["positions"][:]

        dump.write("ITEM: TIMESTEP\n")
        dump.write(f"{step}\n")

        dump.write("ITEM: BOX BOUNDS xx yy zz\n")
        dump.write(f"{x_lo:g} {x_hi:g}\n")
        dump.write(f"{y_lo:g} {y_hi:g}\n")
        dump.write(f"{z_lo:g} {z_hi:g}\n")

        dump.write("ITEM: NUMBER OF ATOMS\n")
        dump.write(f"{len(points)}\n")

        dump.write("ITEM: ATOMS id x y z chain_id\n")

        for chain_index, (beg, end) in enumerate(chain_ranges):
            chain_id = chain_index + 1  # Must start from 1
            for i in range(beg, end):
                atom_id = i + 1  # Must start from 1
                x, y, z = points[i]
                dump.write(f"{atom_id} {x:g} {y:g} {z:g} {chain_id}\n")


def determine_chain_ranges(beads: pd.DataFrame):
    chain_codes = beads["chain"].astype("category").cat.codes
    cur_start = 0
    cur_chain = chain_codes[0]

    chain_ranges = []
    for i, chain in enumerate(chain_codes):
        if chain != cur_chain:
            chain_ranges.append((cur_start, i))
            cur_start = i
            cur_chain = chain
    chain_ranges.append((cur_start, len(chain_codes)))

    return np.array(chain_ranges)


def determine_particle_types(beads: pd.DataFrame):
    particle_types = []
    for a in beads["A"]:
        if a == 1:
            particle_types.append(PARTICLE_TYPE_A)
        elif a == 0:
            particle_types.append(PARTICLE_TYPE_B)
        else:
            particle_types.append(PARTICLE_TYPE_U)
    return np.array(particle_types)
