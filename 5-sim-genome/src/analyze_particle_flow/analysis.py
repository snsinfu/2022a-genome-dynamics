import hashlib
import json
import multiprocessing
import os.path

import h5py
import numba
import numpy as np
import scipy.spatial

from .utils import gaussian_smooth


ANALYSIS_ROOT = "/particle_flow"
HASHNAME_LENGTH = 7


def run(outfile, trajfiles, *, name, smoothing, velocity_delay, scan_radius, jobs):

    config = {
        "smoothing": smoothing,
        "velocity_delay": velocity_delay,
        "scan_radius": scan_radius,
    }
    config_json = json.dumps(config)

    if name is None:
        name = hashlib.sha256(config_json.encode("utf-8")).hexdigest()
        name = name[:HASHNAME_LENGTH]

    with h5py.File(outfile, "a") as output:
        analysis = put_group(output, f"{ANALYSIS_ROOT}/{name}")
        put_dataset(analysis, ".config", config_json)

        sample_names = []

        for trajfile in trajfiles:
            sample_name, _ = os.path.splitext(os.path.basename(trajfile))
            sample_names.append(sample_name)

            with h5py.File(trajfile, "r") as input:
                positions_history = load_positions(input)
            if smoothing:
                positions_history = gaussian_smooth(positions_history, smoothing)
            velocities_history = compute_velocities(positions_history, velocity_delay)
            frame_count = len(positions_history)

            if jobs is None:
                flows_history = []
                for frame in range(frame_count):
                    positions = positions_history[frame]
                    velocities = velocities_history[frame]
                    flows = compute_flow(positions, velocities, scan_radius)
                    flows_history.append(flows)
            else:
                with multiprocessing.Pool(jobs) as pool:
                    flows_history = pool.starmap(
                        compute_flow,
                        zip(
                            positions_history,
                            velocities_history,
                            [scan_radius] * frame_count,
                        ),
                    )

            flows_history = np.array(flows_history, dtype=np.float32)

            put_dataset(
                analysis,
                f"{sample_name}/position",
                positions_history,
                shuffle=True,
                compression=1,
            )

            put_dataset(
                analysis,
                f"{sample_name}/velocity",
                flows_history,
                shuffle=True,
                compression=1,
            )

        # h5py requires variable-length strings to be encoded as bytes.
        put_dataset(
            analysis, ".samples", [name.encode("utf-8") for name in sample_names]
        )


def put_group(store, path):
    if path in store:
        return store[path]
    return store.create_group(path)


def put_dataset(store, path, data=None, **kwargs):
    if path in store:
        del store[path]
    return store.create_dataset(path, data=data, **kwargs)


def compute_flow(points, velocities, scan_radius):
    kdtree = scipy.spatial.cKDTree(points)
    pairs = kdtree.query_pairs(scan_radius, output_type="ndarray")
    del kdtree

    @numba.njit
    def collect(sums, counts, pairs, peer_values):
        for pair_index in range(len(pairs)):
            i, j = pairs[pair_index]
            sums[i] += peer_values[j]
            sums[j] += peer_values[i]
            counts[i] += 1
            counts[j] += 1

    flow_velocities = np.zeros(velocities.shape, dtype=np.float32)
    contacts = np.zeros(len(points), dtype=np.float32)
    collect(flow_velocities, contacts, pairs, velocities)

    # Include the velocity of the central particle.
    flow_velocities += velocities
    contacts += 1

    flow_velocities /= contacts[:, None]

    del pairs
    del contacts

    return flow_velocities


def compute_velocities(positions_history, delay):
    return np.array([
        estimate_velocity(positions_history, frame, delay)
        for frame in range(len(positions_history))
    ])


def estimate_velocity(positions_history, frame, delay):
    """
    Estimate velocity by linear regression on path.
    """
    window = delay + 1
    back_window = window // 2
    forw_window = window - back_window

    if frame - back_window < 0:
        back_window = frame
    if frame + forw_window > len(positions_history):
        forw_window = len(positions_history) - frame

    window = back_window + forw_window

    times = np.arange(window)
    centered_times = times - times.mean(0)
    cov_to_slope = 1 / np.sum(np.square(centered_times))

    paths = positions_history[frame - back_window:frame + forw_window]
    centered_paths = paths - paths.mean(0)

    velocities = np.einsum("t,tik->ik", centered_times, centered_paths) * cov_to_slope

    return velocities


def load_positions(input):
    positions_history = []
    samples = input["snapshots/interphase"]
    for step in samples[".steps"]:
        positions = samples[step]["positions"][:]
        positions_history.append(positions)
    return np.array(positions_history)
