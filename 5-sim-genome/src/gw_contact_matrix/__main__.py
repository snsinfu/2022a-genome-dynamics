import argparse
import os
import signal

from .command import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--frame-range",
        type=str,
        default=None,
        help="specify frames to analyze in start[:end] format",
    )

    parser.add_argument(
        "--rebin-rate",
        type=int,
        default=1,
        help="specify binning rate (default: 1)",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="specify output filename",
    )

    parser.add_argument(
        "inputs",
        type=str,
        nargs="+",
        help="specify trajectory files",
    )

    return cook_args(vars(parser.parse_args()))


def cook_args(args):
    frame_range = args["frame_range"]
    if frame_range is not None:
        tokens = frame_range.split(":")
        if len(tokens) == 1:
            frame_range = int(tokens[0]), None
        if len(tokens) == 2:
            frame_range = int(tokens[0]), int(tokens[1])
        args["frame_range"] = frame_range
    return args


try:
    main()
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
except BrokenPipeError:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
