import argparse
import os
import signal

from .nad_profile import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--after", type=int, default=None)
    parser.add_argument("--before", type=int, default=None)
    parser.add_argument("--chroms", type=str)
    parser.add_argument("jobdir", type=str)
    return postprocess_args(vars(parser.parse_args()))


def postprocess_args(args):
    args["chroms"] = args["chroms"].split(",")
    return args


def untrap(sig):
    signal.signal(sig, signal.SIG_DFL)
    os.kill(os.getpid(), sig)


try:
    main()
except KeyboardInterrupt:
    untrap(signal.SIGINT)
except BrokenPipeError:
    untrap(signal.SIGPIPE)
