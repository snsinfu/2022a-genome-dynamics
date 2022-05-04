import h5py
import numpy as np


def load_contact_matrices(datasets, norm_method="RAW"):
    contact_matrices = {}
    chrom_sizes = {}

    # Bin definitions. We do not need start/end coordinates here.
    chroms = datasets["bins/chrom"][:]
    ends = datasets["bins/end"][:]

    # Normalization vector.
    if norm_method == "RAW":
        norm_vector = np.ones(chroms.shape[0])
    else:
        norm_vector = datasets[f"bins/{norm_method}"][:]

    # Mapping between chromosome name and values in the `chroms` array.
    chrom_name2code = h5py.h5t.check_dtype(enum=chroms.dtype)
    sorted_keys = [
        name for name, code in sorted(chrom_name2code.items(), key=lambda kv: kv[1])
    ]

    # Identify chromosome segments.
    chrom_ranges = {}
    indices = np.arange(chroms.shape[0])

    for key in chrom_name2code:
        chrom_range = indices[chroms == chrom_name2code[key]]
        chrom_ranges[key] = (chrom_range[0], chrom_range[-1] + 1)
        chrom_sizes[key] = ends[chrom_range[-1]]

    # Initialize intra-chromosomal contact matrices.
    for key in sorted_keys:
        start, end = chrom_ranges[key]
        size = end - start
        contact_matrices[key] = np.zeros((size, size), dtype=np.float32)

    # (i,j,v) contact map. This is huge. Do not load this into memory!
    bin1 = datasets["pixels/bin1_id"]
    bin2 = datasets["pixels/bin2_id"]
    count = datasets["pixels/count"]
    assert bin1.shape == bin2.shape == count.shape

    # Extract cis contact matrices.
    n_samples = bin1.shape[0]
    chunk_size = 500_000

    for start_idx in range(0, n_samples, chunk_size):
        end_idx = min(start_idx + chunk_size, n_samples)

        ii = bin1[start_idx:end_idx]
        jj = bin2[start_idx:end_idx]
        cc = count[start_idx:end_idx]
        ii, jj = np.minimum(ii, jj), np.maximum(ii, jj)

        cis_selector = (chroms[ii] == chroms[jj])
        ii = ii[cis_selector]
        jj = jj[cis_selector]
        cc = cc[cis_selector] / (norm_vector[ii] * norm_vector[jj])

        # Here we are looping because the chunked data may cross chromosome
        # boundaries.
        for key in chrom_name2code:
            start, end = chrom_ranges[key]
            chrom_selector = (ii >= start) & (ii < end)

            # Convert global index to intra-chromosomal index.
            local_ii = ii[chrom_selector] - start
            local_jj = jj[chrom_selector] - start
            local_cc = cc[chrom_selector]

            contact_matrix = contact_matrices[key]
            contact_matrix[local_ii, local_jj] += local_cc
            contact_matrix[local_jj, local_ii] += local_cc

    return contact_matrices, chrom_sizes


def compute_enrichment_matrices(contact_matrices, valid_chroms=None):
    mean_contact_profile = _compute_mean_contact_profile(contact_matrices, valid_chroms)
    enrichment_matrices = _compute_enrichment_matrices(contact_matrices, mean_contact_profile)
    return enrichment_matrices, mean_contact_profile


def _compute_mean_contact_profile(contact_matrices, valid_chroms=None):
    if valid_chroms is None:
        valid_chroms = contact_matrices.keys()

    max_size = max(m.shape[0] for m in contact_matrices.values())
    contacts = np.zeros(max_size)
    counts = np.zeros(max_size)

    for chrom in valid_chroms:
        contact_matrix = contact_matrices[chrom].copy()
        contact_matrix[contact_matrix == 0] = np.nan

        for n in range(contact_matrix.shape[0]):
            diag = np.diag(contact_matrix, k=n)
            contacts[n] += np.nansum(diag)
            counts[n] += np.sum(1 - np.isnan(diag))

    return contacts / counts


def _compute_enrichment_matrices(contact_matrices, mean_contact_profile):
    enrichment_matrices = {}

    for chrom in contact_matrices.keys():
        contact_matrix = contact_matrices[chrom]
        indices = np.arange(contact_matrix.shape[0])
        distance_matrix = np.abs(indices[:, None] - indices[None, :])
        enrichment_matrix = contact_matrix / mean_contact_profile[distance_matrix]
        enrichment_matrices[chrom] = enrichment_matrix

    return enrichment_matrices


def compute_contact_pca(matrix, valid_mask=None):
    if valid_mask is None:
        # Can't handle zero-valued row.
        valid_mask = np.any(matrix != 0, axis=1)
    valid_matrix = matrix[valid_mask, :][:, valid_mask]

    center = np.mean(valid_matrix, axis=0)
    centered_matrix = valid_matrix - center[None, :]

    u, s, vh = np.linalg.svd(centered_matrix)
    pcs = u * np.sqrt(centered_matrix.shape[0] - 1)
    axes = vh

    # Fill in masked rows with nans.
    extended_pcs = np.empty((valid_matrix.shape[0], matrix.shape[1]))
    extended_pcs[:, :] = np.nan
    for i in range(extended_pcs.shape[0]):
        extended_pcs[i, valid_mask] = pcs[:, i]
    extended_pcs = extended_pcs.T

    explained_variances = s ** 2

    extended_axes = np.empty((valid_matrix.shape[0], matrix.shape[1]))
    extended_axes[:, :] = np.nan
    for i in range(extended_axes.shape[0]):
        extended_axes[i, valid_mask] = axes[i, :]

    return extended_pcs, explained_variances, extended_axes
