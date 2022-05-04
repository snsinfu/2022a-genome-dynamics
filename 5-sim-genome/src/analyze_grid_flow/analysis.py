import hashlib
import json
import multiprocessing
import os.path

import h5py
import numba
import numpy as np
import scipy.spatial

from .utils import gaussian_smooth, estimate_velocity


ANALYSIS_ROOT = "/grid_flow"
HASHNAME_LENGTH = 7


def run(
    outfile,
    trajfiles,
    *,
    name=None,
    smoothing=None,
    velocity_delay=1,
    scan_radius,
    grid_interval,
    x_range,
    y_range,
    z_range,
    jobs=None,
):
    # We store all parameters that may change the numerical result in the
    # output file, for reproducibility.
    config = {
        "smoothing": smoothing,
        "velocity_delay": velocity_delay,
        "scan_radius": scan_radius,
        "grid_interval": grid_interval,
        "x_range": x_range,
        "y_range": y_range,
        "z_range": z_range,
    }
    config_json = json.dumps(config)

    # Use truncated hash of the config as a unique analysis name if not
    # specified. This way multiple analysises with different parameters
    # can be stored in a single HDF5 file.
    if name is None:
        name = hashlib.sha256(config_json.encode("utf-8")).hexdigest()
        name = name[:HASHNAME_LENGTH]

    # arange is inclusive-exclusive while we want inclusive ends.
    epsilon = grid_interval * 0.1
    x_mesh = np.arange(x_range[0], x_range[1] + epsilon, grid_interval)
    y_mesh = np.arange(y_range[0], y_range[1] + epsilon, grid_interval)
    z_mesh = np.arange(z_range[0], z_range[1] + epsilon, grid_interval)
    grid_points, grid_indices = make_grid(x_mesh, y_mesh, z_mesh)
    grid_kdtree = scipy.spatial.cKDTree(grid_points)
    grid_shape = [len(x_mesh), len(y_mesh), len(z_mesh)]

    with h5py.File(outfile, "a") as output:
        analysis = put_group(output, f"{ANALYSIS_ROOT}/{name}")
        put_dataset(analysis, ".config", config_json)

        # There may already be some samples in the output file in case of
        # incremental analysis.
        if ".samples" in analysis:
            sample_names = [name.decode("utf-8") for name in analysis[".samples"]]
        else:
            sample_names = []

        put_dataset(analysis, ".grid/shape", grid_shape)
        put_dataset(analysis, ".grid/points", grid_points)
        put_dataset(analysis, ".grid/indices", grid_indices)

        for trajfile in trajfiles:
            # We expect unique filename for each sample.
            sample_name, _ = os.path.splitext(os.path.basename(trajfile))
            sample_names.append(sample_name)

            with h5py.File(trajfile, "r") as input:
                positions_history = load_positions(input)
            if smoothing:
                positions_history = gaussian_smooth(positions_history, smoothing)
            velocities_history = compute_velocities(positions_history, velocity_delay)
            frame_count = len(positions_history)

            if jobs is None:
                results = []
                for frame in range(frame_count):
                    positions = positions_history[frame]
                    velocities = velocities_history[frame]
                    results.append(
                        compute_flow(grid_kdtree, positions, velocities, scan_radius)
                    )
            else:
                with multiprocessing.Pool(jobs) as pool:
                    results = pool.starmap(
                        compute_flow,
                        zip(
                            [grid_kdtree] * frame_count,
                            positions_history,
                            velocities_history,
                            [scan_radius] * frame_count,
                        ),
                    )

            flows_history = np.array(
                [flows for flows, _ in results],
                dtype=np.float32,
            )
            coverages_history = np.array(
                [coverages for _, coverages in results],
                dtype=np.int32,
            )
            del results

            put_dataset(
                analysis,
                f"{sample_name}/flows",
                flows_history,
                shuffle=True,
                scaleoffset=estimate_scaleoffset_factor(flows_history, q=1),
                compression=1,
            )

            put_dataset(
                analysis,
                f"{sample_name}/coverages",
                coverages_history,
                shuffle=True,
                scaleoffset=0, # Lossless yet effective for integral dataset
                compression=1,
            )

            output.flush()

        # Duplicates happen when overwriting previous analysis.
        sample_names = remove_duplicates(sample_names)

        # h5py requires variable-length strings to be encoded as bytes.
        put_dataset(
            analysis, ".samples", [name.encode("utf-8") for name in sample_names]
        )


def remove_duplicates(xs):
    """
    Remove duplicate elements in given list. This function preserves the
    order. For duplicates the last one is preserved.
    """
    result = []
    for i, x in enumerate(xs):
        if x not in xs[i + 1:]:
            result.append(x)
    return result


def estimate_scaleoffset_factor(values, q):
    # We truncate digits so that all data above q-percentile must have at
    # least two significant digits preserved.
    resolution = 0.1 * np.percentile(np.abs(values[values > 0]), q=q)
    return -int(np.floor(np.log10(resolution)))


def put_group(store, path):
    if path in store:
        return store[path]
    return store.create_group(path)


def put_dataset(store, path, data=None, **kwargs):
    if path in store:
        del store[path]
    return store.create_dataset(path, data=data, **kwargs)


def compute_flow(test_kdtree, positions, velocities, scan_radius):
    """
    Compute local mean flow around test points.
    """
    test_points = test_kdtree.data
    dimension = 3

    kdtree = scipy.spatial.cKDTree(positions)
    pairs = test_kdtree.sparse_distance_matrix(
        kdtree, scan_radius, output_type="ndarray"
    )
    pairs = np.transpose([pairs["i"], pairs["j"]])
    del kdtree

    @numba.njit
    def collect(sums, counts, pairs, velocities):
        for pair_index in range(len(pairs)):
            i, j = pairs[pair_index]
            sums[i] += velocities[j]
            counts[i] += 1

    flows = np.zeros((len(test_points), dimension), dtype=np.float32)
    coverages = np.zeros(len(test_points), dtype=np.int32)
    collect(flows, coverages, pairs, velocities)
    del pairs

    # Avoid 0/0 = NaN. We store coverage in the output file, so we can always
    # detect zero-sample grid points without abusing NaNs.
    flows /= np.maximum(coverages, 1)[:, None]

    return flows, coverages


def compute_velocities(positions_history, delay):
    return np.array(
        [
            estimate_velocity(positions_history, frame, delay)
            for frame in range(len(positions_history))
        ]
    )


def load_positions(input):
    positions_history = []
    samples = input["snapshots/interphase"]
    for step in samples[".steps"]:
        positions = samples[step]["positions"][:]
        positions_history.append(positions)
    return np.array(positions_history)


def make_grid(x, y, z):
    """
    Create a three-dimensional grid on given coordinates. Return a tuple of
    arrays of grid points and grid indices. The arrays are flattened and thus
    in the shape of (N, 3) where N is the number of grid points.
    """
    dimension = 3
    ix = np.arange(len(x))
    iy = np.arange(len(y))
    iz = np.arange(len(z))

    grid_points = np.moveaxis(np.meshgrid(x, y, z), 0, -1)
    grid_indices = np.moveaxis(np.meshgrid(ix, iy, iz), 0, -1)

    grid_points = grid_points.reshape(-1, dimension)
    grid_indices = grid_indices.reshape(-1, dimension)

    return grid_points, grid_indices
