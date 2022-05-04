import sys

import h5py
import numba
import numpy as np


def run(*, mcoolfile: str, width: int, binsize: int, outfile: str = None):
    if outfile is None:
        output = sys.stdout
    else:
        output = open(outfile, "w")

    with h5py.File(mcoolfile, "r") as mcool:
        dataset = mcool[f"resolutions/{binsize}"]
        contact_band, chrom_codes, chrom_coords = load_contact_band(
            dataset, band_width=width
        )

    chrom_ranges = enumerate_runs(chrom_codes)
    w_sym = compute_W(contact_band, chrom_ranges)

    log_s = np.log(np.arange(1, w_sym.shape[1]))
    log_w = np.log(w_sym[:, 1:])
    alphas = -estimate_slope(log_s, log_w)

    output.write("chrom\tstart\tend\talpha\n")

    for chrom, key in chrom_codes.dtype.metadata["enum"].items():
        start, end = chrom_ranges[key]
        for i in range(start, end):
            start_coord, end_coord = chrom_coords[i]
            output.write(f"{chrom}\t{start_coord}\t{end_coord}\t{alphas[i]:g}\n")


def compute_W(contact_band, chrom_ranges):
    """
    Given contact matrix C(i, j), passed as contact_band, this function
    constructs local contact decay profile

        W(i, s) = (W+(i, s) + W-(i, s)) / 2

    where

        W+(i, s) = C(i, i+s) / sqrt(C(i, i) C(i+s, i+s)) and
        W-(i, s) = C(i, i-s) / sqrt(C(i, i) C(i-s, i-s)) .
    """
    w_forw = np.empty(contact_band.shape, dtype=np.float32)
    w_back = np.empty(contact_band.shape, dtype=np.float32)
    w_forw[:, :] = np.nan
    w_back[:, :] = np.nan

    for start, end in chrom_ranges:
        for s in range(contact_band.shape[1]):
            compute_forward_W(contact_band[start:end], s, w_forw[start:end, s])
            compute_backward_W(contact_band[start:end], s, w_back[start:end, s])

    return np.nanmean([w_forw, w_back], axis=0)


@numba.jit(error_model="numpy")
def compute_forward_W(band, s, W):
    n = len(band)
    for i in range(0, n - s):
        W[i] = band[i, s] / (band[i, 0] * band[i + s, 0]) ** 0.5


@numba.jit(error_model="numpy")
def compute_backward_W(band, s, W):
    n = len(band)
    for i in range(s, n):
        W[i] = band[i - s, s] / (band[i, 0] * band[i - s, 0]) ** 0.5


def load_contact_band(dataset, band_width: int, norm_method: str = None):
    bins = dataset["bins"]
    chrom_codes = bins["chrom"][:]
    chrom_coords = np.transpose([bins["start"][:], bins["end"][:]])
    n_bins = chrom_codes.shape[0]

    if norm_method is not None:
        norm_vector = dataset[f"bins/{norm_method}"][:]

    # (i, s) element stores the (i, i+s) element of the contact matrix.
    contact_band = np.zeros((n_bins, band_width + 1), dtype=np.float32)

    samples = dataset["pixels"]
    samples_bin1_dataset = samples["bin1_id"]
    samples_bin2_dataset = samples["bin2_id"]
    samples_count_dataset = samples["count"]
    n_samples = samples_bin1_dataset.shape[0]
    chunk_size = 500_000

    for start in range(0, n_samples, chunk_size):
        end = min(start + chunk_size, n_samples)

        ii = samples_bin1_dataset[start:end]
        jj = samples_bin2_dataset[start:end]
        cc = samples_count_dataset[start:end]

        # Select intra-chromosome contact samples
        selection = chrom_codes[ii] == chrom_codes[jj]
        ii = ii[selection]
        jj = jj[selection]
        cc = cc[selection]

        # Compute separation
        ss = np.abs(ii - jj)

        # Select diagonal band
        selection = ss <= band_width
        ii = ii[selection]
        jj = jj[selection]
        cc = cc[selection]
        ss = ss[selection]
        ii, jj = np.minimum(ii, jj), np.maximum(ii, jj)

        if norm_method is not None:
            cc = cc.astype(np.float32)
            cc /= norm_vector[ii]
            cc /= norm_vector[jj]

        contact_band[ii, ss] = cc

    # Hi-C should detect nonzero contact near diagonal. So, zero means an
    # unmappable region. Fill these regions with NaNs.
    contact_band[contact_band == 0] = np.nan

    return contact_band, chrom_codes, chrom_coords


def estimate_slope(xs, ys, axis=-1):
    """
    Compute least-squares slope of each of the rows of `xs` and `ys`.
    """
    mx = np.nanmean(xs, axis=axis)
    my = np.nanmean(ys, axis=axis)
    mxx = np.nanmean(xs * xs, axis=axis)
    mxy = np.nanmean(xs * ys, axis=axis)
    return (mxy - mx * my) / (mxx - mx * mx)


def enumerate_runs(xs):
    """
    Return a list of ranges of equal values in `xs`.
    """
    ranges = []
    start = 0
    while start < len(xs):
        end = find_first_not(xs, xs[start], start=(start + 1))
        ranges.append((start, end))
        start = end
    return ranges


def find_first_not(xs, a, start=0):
    """
    Find the first element of `xs` that does not equal to `a`. Return the index
    of the element if exists, or `len(xs)` otherwise.
    """
    for i in range(start, len(xs)):
        if xs[i] != a:
            return i
    return len(xs)
