#
# Output format:
#
# chain   start         end         A          B          tags
# <chain> <coord_start> <coord_end> <A_weight> <B_weight> <tags>
#
# tags:
#   A    - A-type region
#   u    - u-type region
#   B    - B-type region
#   L1   - LINE L1 region
#   cen  - Centromere
#   anor - Active rDNA
#   bnor - Silenced rDNA
#

import argparse
import math
import os
import signal
import sys

import numpy as np
import pandas as pd


# Ordered name of standard autosomes.
STANDARD_AUTOSOMES = [f"chr{n}" for n in range(1, 23)]

# Haploid chromosome chains bearing active NOR.
STANDARD_ACTIVE_NORS = {"chr13:a", "chr14:a", "chr15:a", "chr21:a", "chr22:a"}

# Parameters for chromatin type classification.
AB_WINDOW = 10
AB_CRITERION = "D1"
AB_CRITERION_SIGN = -1

AB_CRITERION_SIGNS = {
    "D1": -1,
    "I1": -1,
}

# Chromatin types.
TYPE_A = 1
TYPE_B = 2
TYPE_U = 3

TYPE_HEURISTICS = [
    ("cen", TYPE_B),
    ("anor", TYPE_A),
    ("bnor", TYPE_B),
    ("het", TYPE_B),
]

TYPE_TAGS = {
    TYPE_A: "A",
    TYPE_B: "B",
    TYPE_U: "u",
}

TYPE_AB_PARAMETERS = {
    TYPE_A: (1.0, 0.0),
    TYPE_B: (0.0, 1.0),
    TYPE_U: (0.5, 0.5),
}

# Parameters for LINE-rich region detection.
LINE1_WINDOW = 10
LINE1_THRESHOLD_QUANTILE = 0.75

# Column names used for validating input files.
GENOMIC_ANNOTATION_COLUMNS = ["chrom", "start", "end", "GC", "L1", "tags"]
CHROMATIN_SIGNAL_COLUMNS = ["chrom", "start", "end", "D1", "D2", "D3", "I1", "I2"]


def run(*, annotfile, signalfile, outfile, **kwargs):
    if outfile is None:
        output_file = sys.stdout
    else:
        output_file = open(outfile, "w")

    annot_table = pd.read_csv(annotfile, sep="\t")
    assert list(annot_table.columns) == GENOMIC_ANNOTATION_COLUMNS

    signal_table = pd.read_csv(signalfile, sep="\t")
    assert list(signal_table.columns) == CHROMATIN_SIGNAL_COLUMNS

    output_file.write("chain\tstart\tend\tA\tB\ttags\n")

    beads = make_quantized_beads(annot_table, signal_table, **kwargs)
    for chain, start, end, a, b, tags in beads:
        tags = ",".join(tags)
        output_file.write(f"{chain}\t{start}\t{end}\t{a:g}\t{b:g}\t{tags}\n")


def make_quantized_beads(
    annot_table,
    signal_table,
    *,
    use_diploid,
    use_xa,
    use_xi,
    use_y,
    criterion,
    threshold,
    hyperactive_nor,
):
    # Transform interaction signal to z-score for chromatin classification.
    score_table = pd.DataFrame(
        {
            "chrom": signal_table.chrom,
            "score": AB_CRITERION_SIGNS[criterion]
            * robust_standardize(
                signal_table[criterion],
                mask=signal_table.chrom.isin(STANDARD_AUTOSOMES),
            ),
        }
    )

    # We annotate LINEL1-rich chromatin for junk-removal simulation.
    annot_table.L1.replace(0, np.nan)
    smoothout_column(annot_table, "L1", window=LINE1_WINDOW)
    line1_threshold = annot_table.L1.quantile(LINE1_THRESHOLD_QUANTILE)

    def determine_chromatin_class(score, tags):
        """
        Determine chromatin class based on score and tags.
        """
        if math.isnan(score):
            return infer_chromatin_class(tags)

        lo = -threshold
        hi = threshold

        if score < lo:
            return TYPE_A
        if score > hi:
            return TYPE_B

        return TYPE_U

    def infer_chromatin_class(tags):
        """
        Heuristically determine the chromatin class from tags.
        """
        for key, typ in TYPE_HEURISTICS:
            if key in tags:
                return typ
        return TYPE_U

    def determine_parameters(chain_name, annot, score):
        tags = list(set(annot.tags.split(",")))

        # Annotate L1-rich region.
        if annot.L1 > line1_threshold:
            tags.append("L1")

        # Annotate active NOR as anor.
        if "nor" in tags:
            tags.remove("nor")
            if hyperactive_nor or chain_name in STANDARD_ACTIVE_NORS:
                tags.append("anor")
            else:
                tags.append("bnor")

        if chain_name == "chrX:b":
            # Special case: Secondary chrX is silenced and thus completely B-type.
            typ = TYPE_B
        else:
            typ = determine_chromatin_class(score, tags)

        # Tags contain "het" and "eu" inherited from cytogenic annotation. We
        # don't need them.
        if "het" in tags:
            tags.remove("het")
        if "eu" in tags:
            tags.remove("eu")

        tags.insert(0, TYPE_TAGS[typ])

        # Score.
        a, b = TYPE_AB_PARAMETERS[typ]

        return a, b, tags

    def define_beads(chrom, flavor, derive_xa=False):
        chain_name = f"{chrom}:{flavor}"
        annot_track = annot_table[annot_table.chrom == chrom]
        score_track = score_table[score_table.chrom == chrom]["score"]

        if derive_xa:
            score_track = robust_standardize(score_track)

        for i in range(annot_track.shape[0]):
            score = score_track.iloc[i]
            annot = annot_track.iloc[i]
            a, b, tags = determine_parameters(chain_name, annot, score)
            yield chain_name, annot.start, annot.end, a, b, tags

    # Diploid chains:
    #   chr1:a, chr2:a, ..., chr22:a, chrX:a,
    #   chr1:b, chr2:b, ..., chr22:b, chrX:b/chrY:b

    for chrom in STANDARD_AUTOSOMES:
        yield from define_beads(chrom, flavor="a")

    if use_xa:
        yield from define_beads("chrX", flavor="a", derive_xa=(use_xa and use_xi))

    if use_diploid:
        for chrom in STANDARD_AUTOSOMES:
            yield from define_beads(chrom, flavor="b")

    if use_xi:
        yield from define_beads("chrX", flavor="b")

    if use_y:
        yield from define_beads("chrY", flavor="b")


def smoothout_column(table, name, window):
    """
    Apply rolling mean to a column in-place, respecting chromosome boundaries.
    """
    smooth_tracks = []
    for _, track in table[name].groupby(table.chrom, sort=False):
        smooth_tracks.append(track.rolling(window, center=True, min_periods=1).mean())
    table[name] = np.concatenate(smooth_tracks)


def robust_standardize(data, mask=None):
    """
    Apply robust standardization to 1D dataset.
    """
    MAD_FACTOR = 1.4826

    if mask is None:
        train_data = data
    else:
        train_data = data[mask]
    robust_center = np.nanmedian(train_data)
    robust_deviation = np.nanmedian(np.abs(train_data - robust_center)) * MAD_FACTOR

    return (data - robust_center) / robust_deviation
