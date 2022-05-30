import json

from collections import OrderedDict

import h5py
import numpy as np
import pandas as pd

from .defaults import DEFAULT_CONFIG
from .system_definition import TYPE_CEN, TYPE_ANOR, TYPE_NUC, TYPE_CODEBOOK
from .system_definition import make_system_definition


SEED_MAX = 999999

DTYPE_INDEX = np.int32
DTYPE_FLOAT = np.float32
DTYPE_TYPE_CODE = np.int8


def run(*, configfile, genomefile, outputfile, seed=None):
    config, random = load_config(configfile, seed)
    genome = pd.read_csv(genomefile, sep="\t")
    system = make_system_definition(genome, config)

    with h5py.File(outputfile, "w") as store:
        create_hierarchy(store)
        save_config(store, config)
        save_ab_factors(store, system)
        save_particle_types(store, system)
        save_chromosome_ranges(store, system)
        save_centromere_ranges(store, system)
        save_nucleolus_ranges(store, system)
        save_nucleolus_bonds(store, system)


def load_config(filename, seed=None):
    """
    Make simulation configuration based on input JSON file.
    """
    with open(filename, "r") as file:
        config = DEFAULT_CONFIG.copy()
        config.update(json.load(file, object_pairs_hook=OrderedDict))

    # Generate a random seed if not specified. We use a PRNG seeded by this to
    # generate any random number needed.
    if seed is None:
        seed = config.setdefault("seed", np.random.randint(SEED_MAX + 1))
    random = np.random.RandomState(seed)

    # Store the actual seed value to the config for reproducibility. Also we
    # derive seed values used in subsequent simulations here.
    config["seed"] = seed
    config["spindle_seed"] = random.randint(SEED_MAX + 1)
    config["interphase_seed"] = random.randint(SEED_MAX + 1)

    return config, random


def create_hierarchy(store):
    """
    Create HDF5 groups.
    """
    store.create_group("metadata")
    store.create_group("snapshots/spindle")
    store.create_group("snapshots/packing")
    store.create_group("snapshots/relaxation")
    store.create_group("snapshots/interphase")


def save_config(store, config):
    store["metadata/config"] = json.dumps(config)


def save_ab_factors(store, system):
    store["metadata/ab_factors"] = np.array(
        [[part.A, part.B] for part in system.particles], dtype=DTYPE_FLOAT
    )


def save_particle_types(store, system):
    enum_map = {name: code for code, _, name in TYPE_CODEBOOK}
    data = np.array([part.type for part in system.particles], dtype=DTYPE_TYPE_CODE)
    store.create_dataset(
        "metadata/particle_types",
        data=data,
        dtype=h5py.special_dtype(enum=(data.dtype, enum_map)),
    )


def save_chromosome_ranges(store, system):
    store["metadata/chromosome_ranges"] = np.array(
        [[chain.start, chain.end] for chain in system.chromatin_chains],
        dtype=DTYPE_INDEX,
    )
    store["metadata/chromosome_ranges"].attrs["keys"] = json.dumps(
        {chain.name: i for i, chain in enumerate(system.chromatin_chains)}
    )


def save_centromere_ranges(store, system):
    store["metadata/centromere_ranges"] = np.array(
        [[chain.cen_start, chain.cen_end] for chain in system.chromatin_chains],
        dtype=DTYPE_INDEX,
    )
    store["metadata/centromere_ranges"].attrs["keys"] = json.dumps(
        {chain.name: i for i, chain in enumerate(system.chromatin_chains)}
    )


def save_nucleolus_ranges(store, system):
    store["metadata/nucleolus_ranges"] = np.array(
        [[span.start, span.end] for span in system.nucleolus_spans], dtype=DTYPE_INDEX
    )
    store["metadata/nucleolus_ranges"].attrs["keys"] = json.dumps(
        {span.name: i for i, span in enumerate(system.nucleolus_spans)}
    )


def save_nucleolus_bonds(store, system):
    store["metadata/nucleolus_bonds"] = np.array(
        system.nucleolus_bonds, dtype=DTYPE_INDEX
    )
