import json
from dataclasses import dataclass

import numpy as np


PATH_CONFIG = "config"
PATH_CHAIN_RANGES = "chain_ranges"
PATH_LOOPS = "loops_history"
PATH_POSITIONS = "positions_history"


@dataclass
class Trajectory:
    config: dict
    chain_ranges: np.array
    loops: np.array
    positions: np.array


def load_trajectory(store):
    length = store[PATH_POSITIONS].shape[0]

    if PATH_LOOPS in store:
        loops = store[PATH_LOOPS][:]
    else:
        loops = np.empty((length, 0, 3), dtype=int)

    return Trajectory(
        config=json.loads(store[PATH_CONFIG][()]),
        chain_ranges=store[PATH_CHAIN_RANGES][()],
        loops=loops,
        positions=store[PATH_POSITIONS][:],
    )
