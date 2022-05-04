import glob
import json
import os
import sys

import h5py
import numpy as np
import scipy.sparse as sparse


TRAJECTORY_FILENAME = "output-*.h5"


def run(*, jobdir, chroms, before, after):
    paths = glob.glob(os.path.join(jobdir, TRAJECTORY_FILENAME))
    sum_contact_matrix = None

    for path in paths:
        with h5py.File(path, "r") as store:
            types = store["metadata/particle_types"]
            types_enum = h5py.check_dtype(enum=types.dtype)
            types = types[:]

            chrom_ranges = store["metadata/chromosome_ranges"]
            chrom_keys = json.loads(chrom_ranges.attrs["keys"])
            chrom_ranges = chrom_ranges[:]

            for chrom in chroms:
                chrom_index = chrom_keys[chrom]
                beg, end = chrom_ranges[chrom_index]

                contact_matrix = collect_contact_matrix(
                    store["snapshots/interphase"], before=before, after=after, chain=(beg, end)
                )

                if sum_contact_matrix is None:
                    sum_contact_matrix = contact_matrix
                else:
                    sum_contact_matrix += contact_matrix

    np.savetxt(sys.stdout, sum_contact_matrix, fmt="%d", delimiter="\t")


def collect_contact_matrix(phase, before, after, chain):
    beg, end = chain
    chain_size = end - beg
    contact_matrix = np.zeros((chain_size, chain_size), dtype=np.int32)

    sample_keys = phase[".steps"]
    sample_steps = [int(key) for key in sample_keys]

    if before is not None:
        sample_steps = [step for step in sample_steps if step < before]
    if after is not None:
        sample_steps = [step for step in sample_steps if step >= after]

    for step in sample_steps:
        sample = phase[str(step)]

        # Contact map is not saved in every sample frame.
        if "contact_map" not in sample:
            continue

        contact_map_data = sample["contact_map"]

        # Contact map is a very large sparse matrix in the (i,j,v) format. Let's
        # load triplets by chunk to reduce memory pressure.
        if contact_map_data.chunks:
            chunk_rows = contact_map_data.chunks[0]
        else:
            chunk_rows = 1024 * 1024

        for chunk_beg in range(0, contact_map_data.shape[0], chunk_rows):
            chunk_end = min(chunk_beg + chunk_rows, contact_map_data.shape[0])

            i, j, c = contact_map_data[chunk_beg:chunk_end, :].T

            # Select relevant rows.
            selector = (i >= beg) & (i < end)
            selector &= (j >= beg) & (j < end)

            if selector.sum() == 0:
                continue

            i = i[selector] - beg
            j = j[selector] - beg
            c = c[selector]

            contacts = sparse.csc_matrix((c, (i, j)), shape=(chain_size, chain_size))
            contact_matrix += contacts

    contact_matrix = contact_matrix + contact_matrix.T
    np.fill_diagonal(contact_matrix, contact_matrix.max())

    return contact_matrix
