import argparse
import os
import signal

from .dump_lammps import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser(
        prog="dump_lammps",
        description="Create LAMMPS data & dump files from simulation trajectory",
    )
    parser.add_argument(
        "--steps",
        dest="step_range",
        metavar="[<start>:]<end>",
        type=str,
        help="Step range to dump",
    )
    parser.add_argument(
        "trajfile", metavar="<trajectory>", type=str, help="input HDF5 file",
    )
    parser.add_argument(
        "datafile", metavar="<data>", type=str, help="output LAMMPS data file",
    )
    parser.add_argument(
        "dumpfile", metavar="<dump>", type=str, help="output LAMMPS dump file",
    )
    args = parser.parse_args()

    args.step_range = parse_step_range(args.step_range)

    return vars(args)


def parse_step_range(arg):
    """
    Parses step range passed to the `--steps` command-line argument. The
    format is `<start>:<end>`. The start part (`<start>:`) can be omitted
    and defaults to zero. Rreturns tuple `(start, end)`.
    """
    if arg is None:
        return None
    fragments = arg.split(":")
    if len(fragments) == 1:
        return 0, int(arg)
    if len(fragments) == 2:
        return int(fragments[0]), int(fragments[1])
    raise Exception("invalid step range")


try:
    main()
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
except BrokenPipeError:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGPIPE)
