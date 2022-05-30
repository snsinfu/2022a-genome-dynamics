import sys

from dataclasses import dataclass
from typing import Dict, Tuple

import h5py
import numpy as np


MIN_CHUNK_SIZE = 500_000
BLACKLISTED_CHROMS = {"MT"}
NAMED_CHROM_RANK = {"X": 1, "Y": 2, "MT": 3, "M": 3}


def run(*, mcoolfile: str, band_width: int, binsize: int, outfile: str = None):
    with h5py.File(mcoolfile, "r") as mcool:
        datasets = mcool[f"resolutions/{binsize}"]
        band_tracks = extract_forward_bands(datasets, band_width)

    # Output
    if outfile is None:
        output = sys.stdout
    else:
        output = open(outfile, "w")

    header = "chrom\tstart\tend\t"
    header += "\t".join(f"D{k}" for k in range(1, band_width))
    header += "\t"
    header += "\t".join(f"I{k}" for k in range(1, band_width - 1))
    header += "\n"
    output.write(header)

    for chrom, track in band_tracks.items():
        chrom_name = ensure_prefix(chrom, "chr")
        starts = track.starts
        ends = track.ends
        bands = track.forward_bands
        decays = compute_local_decays(bands)
        insulations = compute_insulation_ratios(decays)

        # D0 and I0 are trivial. Start from D1 and I1.
        decays = decays[:, 1:]
        insulations = insulations[:, 1:]

        for start, end, decay, insulation in zip(starts, ends, decays, insulations):
            row = f"{chrom_name}\t{start}\t{end}\t"
            row += "\t".join(f"{value:g}" for value in decay)
            row += "\t"
            row += "\t".join(f"{value:g}" for value in insulation)
            row += "\n"
            output.write(row)


def compute_insulation_ratios(decays: np.ndarray) -> np.ndarray:
    """
    Compute insulation ratios I(i,k):

        I(i,k) = D(i,k) / D(i,k+1) .

    Basically I(i,k) ~ C(i,i+k) / C(i,i+k+1) for small k.
    """
    insulation = np.empty((decays.shape[0], decays.shape[1] - 1))
    for k in range(decays.shape[1] - 1):
        insulation[:, k] = decays[:, k] / decays[:, k + 1]
    return insulation


def compute_local_decays(bands: np.ndarray) -> np.ndarray:
    """
    Compute local decay profile D(i,k):

        D(i,k) = (@(i,k) + @(i,-k)) / 2 ,
        @(i,k) = C(i,i+k) / sqrt(C(i,i) C(i+k,i+k)) .

    Basically D(i,k) ~ C(i,i+k) / C(i,i) for small k.
    """
    decays = np.ones(bands.shape)

    if len(bands) <= 1:
        decays[:] = np.nan
        return decays

    for k in range(1, bands.shape[1]):
        diag_0 = bands[:-k, 0]
        diag_k = bands[k:, 0]
        contact = bands[:-k, k]
        decay = contact / np.sqrt(diag_0 * diag_k)
        assert len(decay) == len(bands) - k

        # Symmetrize decay profile because
        sym_decay = np.concatenate(
            [decay[:k], np.nanmean([decay[k:], decay[:-k]], axis=0), decay[-k:],]
        )
        assert len(sym_decay) == len(bands)

        decays[:, k] = sym_decay

    return decays


@dataclass
class ChromTrack:
    starts: np.ndarray
    ends: np.ndarray
    forward_bands: np.ndarray


def extract_forward_bands(
    datasets: h5py.Group, band_width: int
) -> Dict[str, ChromTrack]:
    # Bin definitions.
    chroms = datasets["bins/chrom"][:]
    starts = datasets["bins/start"][:]
    ends = datasets["bins/end"][:]
    chrom_name2code = h5py.check_dtype(enum=chroms.dtype)

    # Contact map. It is huge! Do not load into memory.
    bin1 = datasets["pixels/bin1_id"]
    bin2 = datasets["pixels/bin2_id"]
    count = datasets["pixels/count"]
    assert bin1.shape == bin2.shape == count.shape

    # Extract diagonal bands. Contact map won't fit in memory so we need to work
    # in chunks.
    n_bins = chroms.shape[0]
    n_samples = bin1.shape[0]
    bands = np.zeros((n_bins, band_width))
    chunk_size = determine_read_chunk(count)

    for start_idx in range(0, n_samples, chunk_size):
        end_idx = min(start_idx + chunk_size, n_samples)

        ii = bin1[start_idx:end_idx]
        jj = bin2[start_idx:end_idx]
        cc = count[start_idx:end_idx]

        # Ensure ii < jj.
        ii, jj = np.minimum(ii, jj), np.maximum(ii, jj)

        # Select in diagonal bands.
        band_selector = (jj - ii < band_width) & (chroms[ii] == chroms[jj])
        ii = ii[band_selector]
        jj = jj[band_selector]
        cc = cc[band_selector]

        # Here we only store forward contacts. Hence, the telomere of q-arm has
        # zero contacts.
        bands[ii, jj - ii] += cc

    # Zero contact in diagonal band most likely means unmappable regions. NaN
    # is appropriate.
    bands[bands == 0] = np.nan

    # Split raw bands data into chromosomes.
    tracks = {}

    chrom_names: List[str] = [
        chrom for chrom in chrom_name2code if chrom not in BLACKLISTED_CHROMS
    ]
    chrom_names.sort(key=by_std_chrom_order)

    for chrom in chrom_names:
        selector = chroms == chrom_name2code[chrom]
        tracks[chrom] = ChromTrack(
            starts=starts[selector], ends=ends[selector], forward_bands=bands[selector],
        )

    return tracks


def by_std_chrom_order(chrom: str) -> Tuple[int, int]:
    # Numeral chromosomes precede. Then, sex and other chromosomes follow. We
    # use tuple to represent this two-level ordering.
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    try:
        n = int(chrom)
        return 0, n
    except:
        pass
    return NAMED_CHROM_RANK[chrom], 0


def determine_read_chunk(dataset: h5py.Dataset) -> int:
    if dataset.chunks is None:
        return MIN_CHUNK_SIZE
    chunk, *_ = dataset.chunks

    # Too small chunk hurts sequential read performance.
    if chunk < MIN_CHUNK_SIZE:
        chunk = (MIN_CHUNK_SIZE + chunk - 1) // chunk * chunk

    return chunk


def ensure_prefix(s: str, prefix: str) -> str:
    if s.startswith(prefix):
        return s
    return prefix + s
