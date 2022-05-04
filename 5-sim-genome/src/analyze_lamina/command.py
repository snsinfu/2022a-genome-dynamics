import json
import os

import h5py
import numpy as np


from .geometry import Ellipsoid


H5_OPTIONS_DISTANCE = {
    "shuffle": True,
    "scaleoffset": 3,
    "compression": "gzip",
    "compression_opts": 1,
}

H5_OPTIONS_CONTACT = {
    "shuffle": True,
    "compression": "gzip",
    "compression_opts": 1,
}

H5_OPTIONS_AVERAGE_CONTACT = {
    "shuffle": True,
    "compression": "gzip",
    "compression_opts": 1,
}


def analyze_distance(*, trajfiles, outfile):
    with h5py.File(outfile, "a") as output:
        # We assume metadata is the same for all trajfiles and transfer the
        # metadata to the output file.
        with h5py.File(trajfiles[0], "r") as store:
            config_json = store["metadata/config"][()]
            particle_types = store["metadata/particle_types"][...]
            ab_factors = store["metadata/ab_factors"][...]
            chromosome_ranges = store["metadata/chromosome_ranges"]
            chromosome_names = [None] * len(chromosome_ranges)
            keys = json.loads(chromosome_ranges.attrs["keys"])
            for name, index in keys.items():
                chromosome_names[index] = name
            chromosome_ranges = chromosome_ranges[...]
            chromosome_names = np.array(chromosome_names, "S")

        config = json.loads(config_json)
        recreate_dataset(output, "metadata/simulation_config", data=config_json)
        recreate_dataset(output, "metadata/particle_types", data=particle_types)
        recreate_dataset(output, "metadata/chromosome_ranges", data=chromosome_ranges)
        recreate_dataset(output, "metadata/chromosome_names", data=chromosome_names)

        for trajfile in trajfiles:
            key, _ = os.path.splitext(os.path.basename(trajfile))

            with h5py.File(trajfile, "r") as store:
                distances_history, scale_history = analyze_distances_history(
                    store["snapshots/interphase"]
                )

            # Record absolute distances from the wall.
            dataset = recreate_dataset(
                output,
                f"distance/{key}",
                shape=distances_history.shape,
                dtype=np.float32,
                **H5_OPTIONS_DISTANCE,
            )
            dataset[...] = distances_history


def analyze_distances_history(snapshots):
    distances_history = []
    scale_history = []

    for step in snapshots[".steps"]:
        sample = snapshots[step]
        context = json.loads(sample["context"][()])
        positions = sample["positions"][:]

        scale = context["bead_scale"]
        wall = Ellipsoid(context["wall_semiaxes"])
        distances = wall.distance_from_surface(positions)

        distances_history.append(distances)
        scale_history.append(scale)

    return np.array(distances_history), np.array(scale_history)


def analyze_contact(*, outfile, name, contact_distance=None):
    analyze_contact_uniform(
        outfile,
        name=name,
        contact_distance=contact_distance,
    )


def analyze_contact_uniform(outfile, name, contact_distance):
    average_contacts_history = None
    num_observations = 0

    with h5py.File(outfile, "a") as output:
        for key in output["distance"].keys():
            distances_history = output["distance"][key][...]
            contacts_history = distances_history < contact_distance

            dataset = recreate_dataset(
                output,
                f"contact/{name}/{key}",
                shape=contacts_history.shape,
                dtype=contacts_history.dtype,
                **H5_OPTIONS_CONTACT,
            )
            dataset[...] = contacts_history

            if average_contacts_history is None:
                average_contacts_history = np.zeros(
                    contacts_history.shape, dtype=np.float32
                )
            average_contacts_history += contacts_history
            num_observations += 1

        average_contacts_history /= num_observations

        dataset = recreate_dataset(
            output,
            f"average_contact/{name}",
            shape=average_contacts_history.shape,
            dtype=average_contacts_history.dtype,
            **H5_OPTIONS_AVERAGE_CONTACT,
        )
        dataset[...] = average_contacts_history


def recreate_dataset(group, name, **kwargs):
    if name in group:
        del group[name]
    return group.create_dataset(name, **kwargs)
