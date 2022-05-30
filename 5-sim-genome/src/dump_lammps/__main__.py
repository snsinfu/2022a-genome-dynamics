import argparse
import os
import signal

from .run import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser(
        prog="dump_lammps",
        description="Create lammps data/dump from simulation trajectory",
    )
    parser.add_argument(
        "--phase",
        dest="phasename",
        metavar="<phase>",
        type=str,
        help="simulation phase to dump (initialization | relaxation | simulation)",
    )
    parser.add_argument(
        "--steps",
        dest="step_range",
        metavar="[<start>:]<end>",
        type=str,
        help="step range to dump",
    )
    parser.add_argument(
        "--apart",
        dest="apart",
        metavar="<factor>",
        type=float,
        help="scale the center of mass of each chain by this factor",
    )
    parser.add_argument(
        "--smooth",
        dest="smooth_window",
        metavar="<window>",
        type=int,
        default=None,
        help="apply gaussian smoothing of given window size",
    )
    parser.add_argument(
        "trajfile", metavar="<trajectory>", type=str, help="input hdf5 file",
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
    if arg is None:
        return None
    fragments = arg.split(":")
    if len(fragments) == 1:
        return 0, int(arg)
    if len(fragments) == 2:
        return int(fragments[0]), int(fragments[1])
    raise Exception("invalid step range")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        os.kill(os.getpid(), signal.SIGINT)
    except BrokenPipeError:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
        os.kill(os.getpid(), signal.SIGPIPE)
