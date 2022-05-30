import argparse
import os
import signal

from .run import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser(
        prog="prepare_trajectory",
        description="Prepare a trajectory file from simulation config and data.",
    )
    parser.add_argument(
        "--seed", metavar="<seed>", type=int, default=None, help="random seed",
    )
    parser.add_argument(
        "configfile", metavar="<config>", type=str, help="simulation config json file",
    )
    parser.add_argument(
        "genomefile", metavar="<genome>", type=str, help="genome definition file",
    )
    parser.add_argument(
        "outputfile", metavar="<output>", type=str, help="output hdf5 file",
    )
    return vars(parser.parse_args())


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        os.kill(os.getpid(), signal.SIGINT)
    except BrokenPipeError:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
        os.kill(os.getpid(), signal.SIGPIPE)
