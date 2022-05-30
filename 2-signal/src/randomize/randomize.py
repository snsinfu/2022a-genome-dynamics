import sys

import numpy as np
import pandas as pd


STRUCTURAL_ELEMENTS = ["cen", "anor", "bnor"]


def run(*, infile, outfile, seed, preserve_structure, completely_random):
    if preserve_structure and completely_random:
        raise Exception(
            "--preserve-structure and --completely-random options can not "
            "both be specified"
        )

    random = np.random.RandomState(seed=seed)
    ab_table = pd.read_csv(infile, sep="\t")

    if outfile is None:
        output = sys.stdout
    else:
        output = open(outfile, "w")

    header = "\t".join(ab_table.columns) + "\n"
    output.write(header)

    order = np.copy(ab_table.index)

    if preserve_structure:
        selector = sum(ab_table.tags.str.contains(x) for x in STRUCTURAL_ELEMENTS) == 0
        # Can't shuffle subarray directly. Shuffle a copy of the indices of
        # non-structural elements and assign the result back.
        sub_order = order[selector]
        random.shuffle(sub_order)
        order[selector] = sub_order
    else:
        random.shuffle(order)

    ab_table.A.values[:] = np.copy(ab_table.A[order])
    ab_table.B.values[:] = np.copy(ab_table.B[order])

    if completely_random:
        ab_table.A.values[:] = random.randint(2, size=len(ab_table))
        ab_table.B.values[:] = 1 - ab_table.A.values

    for i in ab_table.index:
        # Fix A/B tagging
        tags = ab_table.tags[i]
        delta_ab = ab_table.A[i] - ab_table.B[i]

        if delta_ab > 0:
            tags = tags.replace("B", "A")
            tags = tags.replace("u", "A")
        elif delta_ab < 0:
            tags = tags.replace("A", "B")
            tags = tags.replace("u", "B")
        else:
            tags = tags.replace("A", "u")
            tags = tags.replace("B", "u")

        ab_table.tags.values[i] = tags

        print("\t".join(str(x) for x in ab_table.iloc[i]), file=output)
