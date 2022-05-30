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
    sum_contact_profile = None

    for path in paths:
        with h5py.File(path, "r") as store:
            chrom_ranges = store["metadata/chromosome_ranges"]
            chrom_keys = json.loads(chrom_ranges.attrs["keys"])
            chrom_ranges = chrom_ranges[:]

            for chrom in chroms:
                chrom_index = chrom_keys[chrom]
                beg, end = chrom_ranges[chrom_index]

                contact_profile = collect_nucleolus_contacts(
                    store["snapshots/interphase"],
                    before=before, after=after, chain=(beg, end)
                )

                if sum_contact_profile is None:
                    sum_contact_profile = contact_profile
                else:
                    sum_contact_profile += contact_profile

    np.savetxt(sys.stdout, sum_contact_profile, fmt="%d", delimiter="\t")


def collect_nucleolus_contacts(phase, before, after, chain):
    beg, end = chain
    chain_size = end - beg
    contact_profile = np.zeros(chain_size, dtype=np.int32)

    if "metadata" in phase:
        metadata = phase["metadata"]
    else:
        metadata = phase.file["metadata"]

    if "particle_types" not in metadata:
        metadata = phase.file["metadata"]

    types = metadata["particle_types"]
    types_enum = h5py.check_dtype(enum=types.dtype)
    types = types[:]
    nucleolus_type = types_enum["nucleolus"]

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

            # Chromosome vs nucleolus
            selector = (i >= beg) & (i < end) & (types[j] == nucleolus_type)
            contact_profile[i[selector] - beg] += c[selector]

            # Nucleolus vs chromosome
            selector = (j >= beg) & (j < end) & (types[i] == nucleolus_type)
            contact_profile[j[selector] - beg] += c[selector]

    return contact_profile
