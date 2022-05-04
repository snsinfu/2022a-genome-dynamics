import argparse
import os
import signal

from .analysis import run

PROGNAME, _, _ = __spec__.name.rpartition(".")


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser(prog=PROGNAME)
    parser.add_argument("--name", type=str, default=None)
    parser.add_argument("--smoothing", type=int, default=None)
    parser.add_argument("--velocity-delay", type=int, default=1)
    parser.add_argument("--scan-radius", type=float)
    parser.add_argument("--jobs", type=int, default=None)
    parser.add_argument("outfile", type=str)
    parser.add_argument("trajfiles", type=str, nargs="+")
    return vars(parser.parse_args())


def untrap(sig):
    signal.signal(sig, signal.SIG_DFL)
    os.kill(os.getpid(), sig)


try:
    main()
except KeyboardInterrupt:
    untrap(signal.SIGINT)
except BrokenPipeError:
    untrap(signal.SIGPIPE)
