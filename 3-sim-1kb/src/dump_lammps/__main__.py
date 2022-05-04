import argparse
import os
import signal

import h5py

from .lammps import write_data, write_dump
from .trajectory import load_trajectory


def main(*, trajfile, start, end):
    with h5py.File(trajfile, "r") as store:
        basename, _ = os.path.splitext(trajfile)
        datafile = basename + ".data"
        dumpfile = basename + ".dump"

        trajectory = load_trajectory(store)

        if start is None:
            start = 0
        if end is None:
            end = trajectory.positions.shape[0]

        with open(datafile, "w") as data:
            write_data(data, trajectory)

        with open(dumpfile, "w") as dump:
            write_dump(dump, trajectory, start, end)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("trajfile", type=str)
    return vars(parser.parse_args())


try:
    main(**parse_args())
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
