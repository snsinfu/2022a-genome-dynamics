import argparse
import os
import signal

from .model_genome import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        dest="outfile",
        type=str,
        default=None,
        help="output file name (default: stdout)",
    )

    parser.add_argument(
        "--diploid",
        dest="use_diploid",
        action="store_true",
        default=False,
        help="use homologue autosomes",
    )

    parser.add_argument(
        "--Xa",
        dest="use_xa",
        action="store_true",
        default=False,
        help="use activated X chromosome",
    )

    parser.add_argument(
        "--Xi",
        dest="use_xi",
        action="store_true",
        default=False,
        help="use inactivated X chromosome",
    )

    parser.add_argument(
        "--Y",
        dest="use_y",
        action="store_true",
        default=False,
        help="use Y chromosome",
    )

    parser.add_argument(
        "--hyperactive-nor",
        action="store_true",
        default=False,
        help="activate all NORs",
    )

    parser.add_argument(
        "--criterion", default="D1", help="criterion of chromatin classification",
    )

    parser.add_argument(
        "--threshold", type=float, default=0.5, help="A/B threshold z-score",
    )

    parser.add_argument(
        "annotfile", type=str, help="input genome annotation file",
    )

    parser.add_argument(
        "signalfile", type=str, help="input chromatin signal file",
    )

    return vars(parser.parse_args())


try:
    main()
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
except BrokenPipeError:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGPIPE)
