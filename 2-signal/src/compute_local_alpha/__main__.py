import argparse
import os
import signal
import warnings

from .command import run


def main():
    # np.nanmean() warns all-nan input. The function still works and returns
    # nan, which is the expected behavior for us.
    warnings.filterwarnings(
        "ignore", category=RuntimeWarning, message="Mean of empty slice",
    )
    run(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", dest="width", type=int, default=10)
    parser.add_argument("-b", dest="binsize", type=int, default=100000)
    parser.add_argument("-o", dest="outfile", type=str, default=None)
    parser.add_argument("mcoolfile", type=str)
    return vars(parser.parse_args())


try:
    main()
except KeyboardInterrupt:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGINT)
except BrokenPipeError:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGPIPE)
