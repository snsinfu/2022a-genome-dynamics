import json

import h5py
import numpy as np

from .refinement import get_refinement_method


def run(*, trajfile):
    with h5py.File(trajfile, "r+") as store:
        metadata = store["metadata"]

        config = json.loads(metadata["config"][()])
        coarse = config["init_coarse_graining"]
        method = config["init_refinement_method"]
        do_refine = get_refinement_method(method)

        particle_count = metadata["particle_types"].shape[0]
        chrom_ranges = metadata["chromosome_ranges"][:]
        nucleo_bonds = metadata["nucleolus_bonds"][:]

        init = store["snapshots/packing"]
        init_chains = init["metadata/chromosome_ranges"][:]
        init_positions = init[init[".steps"][-1]]["positions"][:]

        fine_positions = np.empty((particle_count, 3))

        # Refine chromosomes.
        for i in range(len(chrom_ranges)):
            beg, end = init_chains[i]
            init_chain = init_positions[beg:end]

            fine_chain = do_refine(init_chain, len(init_chain) * coarse)

            beg, end = chrom_ranges[i]
            fine_positions[beg:end] = fine_chain[:(end - beg)]

        # Attach nucleolar particles.
        for nor, nuc in nucleo_bonds:
            fine_positions[nuc] = fine_positions[nor]

        # Save as the initial conformation of relax phase simulation.
        phase = store["snapshots/relaxation"]
        if "0" in phase:
            del phase["0"]
        phase["0/positions"] = fine_positions
