import json
import sys

import h5py
import numba
import numpy as np
import scipy.sparse


contact_chunk = 1_000_000


def run(*, inputs, output, frame_range, rebin_rate):
    with h5py.File(inputs[0], "r") as store:
        rebin_map, binned_ranges = determine_rebin_map(
            store, rebin_rate
        )

    n_bins = binned_ranges.max()
    contact_matrix = np.zeros((n_bins, n_bins), dtype=np.int32)

    trace("Loading: ")
    progress = load_contact_matrix_into(
        inputs, contact_matrix, rebin_map, frame_range=frame_range
    )
    for i in progress:
        if i % 10 == 0:
            trace(str(i))
        trace(".")
    trace(" DONE\n")

    save_contact_matrix(output, contact_matrix, binned_ranges, rebin_map)


def trace(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()


def save_contact_matrix(filename, contact_matrix, binned_ranges, rebin_map):
    with h5py.File(filename, "w") as store:
        ranges_dataset = store.create_dataset(
            "metadata/chromosome_ranges",
            shape=binned_ranges.shape,
            dtype=binned_ranges.dtype,
        )
        ranges_dataset[...] = binned_ranges

        store["metadata/rebin_map"] = rebin_map

        matrix_dataset = store.create_dataset(
            "contact_matrix",
            shape=contact_matrix.shape,
            dtype=np.int32,
            chunks=True,
            scaleoffset=0,
            shuffle=True,
            compression="gzip",
            compression_opts=1,
        )
        matrix_dataset[...] = contact_matrix


def load_contact_matrix_into(inputs, contact_matrix, rebin_map, frame_range=None):
    if frame_range is None:
        frame_range = None, None

    for input_index, filename in enumerate(inputs):
        with h5py.File(filename, "r") as store:
            snapshots = store["snapshots/interphase"]
            for step in snapshots[".steps"][slice(*frame_range)]:
                sample = snapshots[step]
                if "contact_map" not in sample:
                    continue

                for chunk in iterate_chunks(sample["contact_map"], contact_chunk):
                    collect_contacts(chunk, contact_matrix, rebin_map)
        yield input_index


def iterate_chunks(arr, chunk_size):
    for chunk_start in range(0, len(arr), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(arr))
        yield arr[chunk_start:chunk_end]


@numba.jit
def collect_contacts(samples, matrix, index_map):
    for sample_index in range(len(samples)):
        i, j, v = samples[sample_index]
        n = len(index_map)

        # Ignore non-chromatin contacts.
        if i >= n or j >= n:
            continue

        bin_i = index_map[i]
        bin_j = index_map[j]
        matrix[bin_i, bin_j] += v
        matrix[bin_j, bin_i] += v


def determine_rebin_map(store, rebin_rate):
    src_ranges = store["metadata/chromosome_ranges"]
    chrom_names = json.loads(src_ranges.attrs["keys"])
    src_ranges = src_ranges[:]

    n_src_bins = src_ranges.max()
    rebin_map = np.zeros(n_src_bins, dtype=np.int32)
    binned_ranges = []
    chrom_start = 0

    for start, end in src_ranges:
        bins = np.arange(end - start) // rebin_rate
        rebin_map[start:end] = bins + chrom_start
        chrom_end = chrom_start + bins[-1] + 1
        binned_ranges.append((chrom_start, chrom_end))

        chrom_start = chrom_end
        chrom_end = None

    binned_ranges = np.array(
        binned_ranges,
        dtype=h5py.enum_dtype(chrom_names, np.int32),
    )
    return rebin_map, binned_ranges
