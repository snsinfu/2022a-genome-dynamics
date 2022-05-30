import argparse

from .store import Store


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajfile")
    args = parser.parse_args()

    with Store(args.trajfile, "r") as store:
        interphase = store.phase("interphase")
        for sample in interphase.snapshots:
            print(sample.context)
        print(interphase.metadata)


main()
