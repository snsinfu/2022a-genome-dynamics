import json

import h5py
import numpy as np

from .utils import gaussian_smooth


x_lo, x_hi = -10, 10
y_lo, y_hi = -10, 10
z_lo, z_hi = -10, 10


def run(
    *,
    trajfile,
    datafile,
    dumpfile,
    phasename=None,
    step_range=None,
    apart=None,
    smooth_window=None,
):
    if phasename is None:
        phasename = "simulation"

    with h5py.File(trajfile, "r") as store:
        phase = store["snapshots"][phasename]
        with open(datafile, "w") as data:
            create_data(phase, data)
        with open(dumpfile, "w") as dump:
            create_dump(phase, dump, step_range, apart=apart, smooth_window=smooth_window)


def create_data(phase, data):
    init_step = phase[".steps"][0]
    init_points = phase[init_step]["positions"][:]

    metadata = get_metadata(phase)
    chain_ranges = metadata["chromosome_ranges"][:]

    if "particle_types" in metadata:
        particle_types = metadata["particle_types"][:]
    else:
        particle_types = np.ones(len(init_points))

    n_atoms = len(init_points)
    n_bonds = sum(end - beg - 1 for beg, end in chain_ranges)
    n_atom_types = particle_types.max() - particle_types.min() + 1
    n_bond_types = 1

    data.write("Simulation\n")
    data.write("\n")
    data.write(f"{ n_atoms } atoms\n")
    data.write(f"{ n_bonds } bonds\n")
    data.write("\n")
    data.write(f"{ n_atom_types } atom types\n")
    data.write(f"{ n_bond_types } bond types\n")
    data.write("\n")
    data.write(f"{ x_lo } { x_hi } xlo xhi\n")
    data.write(f"{ y_lo } { x_hi } xlo xhi\n")
    data.write(f"{ z_lo } { x_hi } xlo xhi\n")
    data.write("\n")
    data.write("Atoms\n")
    data.write("\n")

    for i, point in enumerate(init_points):
        atom_id = i + 1
        atom_type = particle_types[i]
        x, y, z = point
        data.write(f"{ atom_id } { atom_type } { x :g} { y :g} { z :g}\n")

    data.write("\n")
    data.write("Bonds\n")
    data.write("\n")

    bond_id = 1
    bond_type = 1
    for beg, end in chain_ranges:
        for i in range(beg, end - 1):
            j = i + 1
            atom_i = i + 1
            atom_j = j + 1
            data.write(f"{ bond_id } { bond_type } { atom_i } { atom_j }\n")
            bond_id += 1


def create_dump(phase, dump, step_range, apart=None, smooth_window=None):
    metadata = get_metadata(phase)
    chain_ranges = metadata["chromosome_ranges"][:]

    points_history = [
        phase[step]["positions"][:] for step in phase[".steps"]
    ]

    if smooth_window:
        points_history = gaussian_smooth(np.array(points_history), window=smooth_window)

    for frame_index, step in enumerate(phase[".steps"]):
        if step_range is not None:
            if int(step) < step_range[0]:
                continue
            if int(step) > step_range[1]:
                continue

        points = points_history[frame_index]

        dump.write("ITEM: TIMESTEP\n")
        dump.write(f"{ step }\n")

        dump.write("ITEM: BOX BOUNDS xx yy zz\n")
        dump.write(f"{ x_lo :g} { x_hi :g}\n")
        dump.write(f"{ y_lo :g} { y_hi :g}\n")
        dump.write(f"{ z_lo :g} { z_hi :g}\n")

        dump.write("ITEM: NUMBER OF ATOMS\n")
        dump.write(f"{ len(points) }\n")

        dump.write("ITEM: ATOMS id x y z chain_id\n")

        for chain_index, (beg, end) in enumerate(chain_ranges):
            chain_id = chain_index + 1

            center = np.zeros(3)
            if apart is not None:
                for i in range(beg, end):
                    center += points[i]
                center /= end - beg
                center *= apart - 1

            for i in range(beg, end):
                atom_id = i + 1
                x, y, z = points[i] + center
                dump.write(f"{ atom_id } { x :g} { y :g} { z :g} { chain_id }\n")

        chain_id += 1
        for i in range(end, len(points)):
            atom_id = i + 1
            x, y, z = points[i]
            dump.write(f"{ atom_id } { x :g} { y :g} { z :g} { chain_id }\n")


def get_metadata(phase):
    if "metadata" in phase:
        return phase["metadata"]
    else:
        return phase.file["metadata"]
