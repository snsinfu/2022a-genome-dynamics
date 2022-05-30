import argparse
import os
import signal

from .command import analyze_distance, analyze_contact


PROGNAME, _, _ = __spec__.name.rpartition(".")


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser(prog=PROGNAME)
    subparsers = parser.add_subparsers(dest="subcommand")

    sub_distance = subparsers.add_parser("distance")
    sub_distance.add_argument("outfile", type=str)
    sub_distance.add_argument("trajfiles", type=str, nargs="+")

    sub_contact = subparsers.add_parser("contact")
    sub_contact.add_argument("--name", type=str, default="uniform")
    sub_contact.add_argument("--contact-distance", type=float, default=None)
    sub_contact.add_argument("outfile", type=str)

    return vars(parser.parse_args())


def run(subcommand, **kwargs):
    if subcommand == "distance":
        analyze_distance(**kwargs)

    if subcommand == "contact":
        analyze_contact(**kwargs)


def untrap(sig):
    signal.signal(sig, signal.SIG_DFL)
    os.kill(os.getpid(), sig)


try:
    main()
except KeyboardInterrupt:
    untrap(signal.SIGINT)
except BrokenPipeError:
    untrap(signal.SIGPIPE)
