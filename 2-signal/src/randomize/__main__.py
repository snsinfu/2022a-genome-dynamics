#
# Script for randomizing AB classification of genome.
#

import argparse
import os
import signal

from . import randomize


def main():
    randomize.run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--preserve-structure", action="store_true", default=False)
    parser.add_argument("--completely-random", action="store_true", default=False)
    parser.add_argument("-o", dest="outfile", type=str, default=None)
    parser.add_argument("infile", type=str)
    return vars(parser.parse_args())


try:
    main()
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
except BrokenPipeError:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGPIPE)
