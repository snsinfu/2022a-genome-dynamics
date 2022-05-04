import argparse
import os
import signal

from .run import run


def main():
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser(
        prog="refine_initialization",
        description="Refine coarse initial conformation",
    )
    parser.add_argument(
        "trajfile", metavar="<trajectory>", type=str, help="in/out trajectory file",
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
