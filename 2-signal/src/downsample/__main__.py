#
# Script for downsampling decay signal data with symmetric boxcar window.
#

import argparse
import os
import signal
import sys

import numpy as np
import pandas as pd


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rate", type=int, default=2)
    parser.add_argument("--window", type=int, default=None)
    parser.add_argument("-o", dest="outfile", type=str, default=None)
    parser.add_argument("infile", type=str)
    return vars(parser.parse_args())


def run(*, infile, outfile, rate, window):
    decay_table = pd.read_csv(infile, sep="\t")

    if outfile is None:
        output = sys.stdout
    else:
        output = open(outfile, "w")

    header = "\t".join(decay_table.columns) + "\n"
    output.write(header)

    for chrom, track in decay_table.groupby(decay_table["chrom"], sort=False):
        starts = track["start"].values
        ends = track["end"].values
        signal = downsample(track.iloc[:, 3:], rate=rate, window=window)

        for i in range(0, len(track), rate):
            start = starts[i]
            end = ends[min(i + rate, len(track)) - 1]
            values = "\t".join(f"{x:g}" for x in signal[i // rate, :])
            row = f"{chrom}\t{start}\t{end}\t{values}\n"
            output.write(row)


def downsample(data, rate=2, window=None):
    """
    Downsample the first axis of 2D array.
    """
    if window is None:
        window = rate
    pad_data = np.pad(
        data, pad_width=[(0, rate), (0, 0)], mode="constant", constant_values=np.nan
    )
    smooth_data = rolling_mean(pad_data, window)
    down_data = smooth_data[rate::rate]
    assert len(down_data) == (len(data) + rate - 1) // rate
    return down_data


def rolling_mean(data, window):
    """
    Computes backward rolling mean of the first axis of 2D array.
    """
    return pd.DataFrame(data).rolling(window, min_periods=1).mean().values


try:
    main()
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
except BrokenPipeError:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGPIPE)
