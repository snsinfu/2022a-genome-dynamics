from dataclasses import dataclass

import numpy as np


@dataclass
class Atom:
    id: int
    type: int
    x: float
    y: float
    z: float
    chain_id: int


@dataclass
class Bond:
    id: int
    type: int
    atom_1: int
    atom_2: int


ATOM_TYPE_MONOMER = 1
ATOM_TYPE_LOOP = 2
BOND_TYPE_CHAIN = 1
BOND_TYPE_LOOP = 2


class LAMMPSCoder:
    def __init__(self, trajectory):
        self._trajectory = trajectory

        # Define topology adjusted for LAMMPS output.
        n_frames, n_monomers, _ = trajectory.positions.shape
        n_loops = trajectory.loops.shape[1]

        self._n_frames = n_frames
        self._n_monomers = n_monomers
        self._n_loops = n_loops
        self._n_loop_factors = n_loops * 2
        self._n_chain_bonds = sum(
            end - start - 1 for start, end in self._trajectory.chain_ranges
        )

        # Starting ID of each type of atoms.
        self._atom_start_monomers = 1
        self._atom_start_loop_factors = 1 + self._n_monomers

        # Starting ID of each type of bonds.
        self._bond_start_chains = 1
        self._bond_start_loops = 1 + self._n_chain_bonds


    def query_atoms(self, frame):
        positions = self._trajectory.positions[frame]
        loops = self._trajectory.loops[frame]

        chain_ids = np.empty(self._n_monomers, dtype=np.int32)
        for i, (start, end) in enumerate(self._trajectory.chain_ranges):
            chain_ids[start:end] = i + 1

        # Chain monomers
        for i in range(self._n_monomers):
            yield Atom(
                id=(self._atom_start_monomers + i),
                type=ATOM_TYPE_MONOMER,
                x=positions[i, 0],
                y=positions[i, 1],
                z=positions[i, 2],
                chain_id=chain_ids[i],
            )

        # Loop factors (two factors per loop)
        for i in range(self._n_loops):
            monomer_1, monomer_2, _ = loops[i]

            if monomer_1 < self._n_monomers and monomer_2 < self._n_monomers:
                pos_1 = positions[monomer_1]
                pos_2 = positions[monomer_2]
                chain_1 = chain_ids[monomer_1]
                chain_2 = chain_ids[monomer_2]
            else:
                pos_1 = [100, 100, 100] # FIXME
                pos_2 = [100, 100, 100] # FIXME
                chain_1 = 0
                chain_2 = 0

            yield Atom(
                id=(self._atom_start_loop_factors + 2 * i),
                type=ATOM_TYPE_LOOP,
                x=pos_1[0],
                y=pos_1[1],
                z=pos_1[2],
                chain_id=chain_1,
            )
            yield Atom(
                id=(self._atom_start_loop_factors + 2 * i + 1),
                type=ATOM_TYPE_LOOP,
                x=pos_2[0],
                y=pos_2[1],
                z=pos_2[2],
                chain_id=chain_2,
            )

    def query_bonds(self):
        # Chain bonds
        bond_id = self._bond_start_chains
        for start, end in self._trajectory.chain_ranges:
            for i in range(start, end - 1):
                yield Bond(
                    id=bond_id,
                    type=BOND_TYPE_CHAIN,
                    atom_1=(self._atom_start_monomers + i),
                    atom_2=(self._atom_start_monomers + i + 1),
                )
            bond_id = bond_id + 1

        # Loop bonds
        for i in range(self._n_loops):
            yield Bond(
                id=(self._bond_start_loops + i),
                type=BOND_TYPE_LOOP,
                atom_1=(self._atom_start_loop_factors + 2 * i),
                atom_2=(self._atom_start_loop_factors + 2 * i + 1),
            )

    def query_bounding_box(self):
        config = self._trajectory.config
        if config["chain"].get("box_size", 0) > 0:
            w = config["chain"]["box_size"] / 2
            min_point = np.array([-w, -w, -w])
            max_point = np.array([w, w, w])
        else:
            dimension = self._trajectory.positions.shape[2]
            points = self._trajectory.positions.reshape(-1, dimension)
            min_point = points.min(axis=0)
            max_point = points.max(axis=0)
        return min_point, max_point

    @property
    def frame_count(self):
        return self._n_frames

    @property
    def atom_count(self):
        return self._n_monomers + self._n_loop_factors

    @property
    def bond_count(self):
        return self._n_chain_bonds + self._n_loops

    @property
    def atom_types(self):
        # (1) Chain monomers and (2) loop factors.
        return 2

    @property
    def bond_types(self):
        # (1) Chain bonds and (2) loops.
        return 2


def write_data(output, trajectory):
    """
    Write LAMMPS data to stream.
    """
    coder = LAMMPSCoder(trajectory)

    atoms = list(coder.query_atoms(0))
    bonds = list(coder.query_bonds())
    min_point, max_point = coder.query_bounding_box()
    min_x, min_y, min_z = min_point
    max_x, max_y, max_z = max_point

    write = lambda s: print(s, file=output)

    write("Simulation data")
    write("")
    write(f"{coder.atom_count} atoms")
    write(f"{coder.bond_count} bonds")
    write("")
    write(f"{coder.atom_types} atom types")
    write(f"{coder.bond_types} bond types")
    write("")
    write(f"{min_x:g} {max_x:g} xlo xhi")
    write(f"{min_y:g} {max_y:g} ylo yhi")
    write(f"{min_z:g} {max_z:g} zlo zhi")

    write("")
    write("Atoms")
    write("")
    for atom in atoms:
        write(f"{atom.id} {atom.type} {atom.x:g} {atom.y:g} {atom.z:g}")

    write("")
    write("Bonds")
    write("")
    for bond in bonds:
        write(f"{bond.id} {bond.type} {bond.atom_1} {bond.atom_2}")


def write_dump(output, trajectory, start, end):
    """
    Write LAMMPS dump to stream.
    """
    coder = LAMMPSCoder(trajectory)

    if start < 0 or start >= coder.frame_count:
        raise Exception("invalid value for the start of frame range")

    if end > coder.frame_count:
        raise Exception("invalid value for the end of frame range")

    min_point, max_point = coder.query_bounding_box()
    min_x, min_y, min_z = min_point
    max_x, max_y, max_z = max_point

    write = lambda s: print(s, file=output)

    for frame in range(start, end):
        write("ITEM: TIMESTEP")
        write(frame)

        write("ITEM: BOX BOUNDS xx yy zz")
        write(f"{min_x:g} {max_x:g}")
        write(f"{min_y:g} {max_y:g}")
        write(f"{min_z:g} {max_z:g}")

        write("ITEM: NUMBER OF ATOMS")
        write(f"{coder.atom_count}")

        write("ITEM: ATOMS id x y z chain_id")
        for atom in coder.query_atoms(frame):
            write(f"{atom.id} {atom.x:g} {atom.y:g} {atom.z:g} {atom.chain_id:d}")
