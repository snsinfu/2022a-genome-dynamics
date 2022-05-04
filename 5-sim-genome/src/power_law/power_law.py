import json

import h5py
import numpy as np

from scipy import sparse
from sklearn.linear_model import LinearRegression


NEAR_RANGE = 3, 20
LONG_RANGE = 20, 100
FAR_RANGE = 100, 1500


def run(trajfiles):
    for trajfile in trajfiles:
        run_once(trajfile)


def run_once(trajfile):
    with h5py.File(trajfile, "r") as store:
        contact_profile = compute_contact_profile(store)
    distances = np.arange(len(contact_profile))

    beg, end = NEAR_RANGE
    near_exp, _ = fit_power_law(distances[beg:end], contact_profile[beg:end])
    beg, end = LONG_RANGE
    long_exp, _ = fit_power_law(distances[beg:end], contact_profile[beg:end])
    beg, end = FAR_RANGE
    far_exp, _ = fit_power_law(distances[beg:end], contact_profile[beg:end])

    print(f"{near_exp:g}\t{long_exp:g}\t{far_exp:g}")


def compute_contact_profile(store):
    particle_types = store["metadata/particle_types"]
    chromosome_ranges = store["metadata/chromosome_ranges"]

    # Exclude non-chromosome particles from analysis by setting NaN
    # chain_ids for these particles.
    chain_ids = np.empty(particle_types.shape[0], dtype=np.float32)
    chain_ids[:] = np.nan

    max_chain_size = 0
    for i, (beg, end) in enumerate(chromosome_ranges):
        chain_ids[beg:end] = i
        max_chain_size = max(max_chain_size, end - beg)

    interphase = store["snapshots/interphase"]
    steps = interphase[".steps"]

    for step in reversed(steps):
        sample = interphase[step]
        if "contact_map" in sample:
            contacts = sample["contact_map"]
            break

    return collect_contact_profile(contacts, chain_ids, size=max_chain_size)


def collect_contact_profile(contacts, chain_ids, *, size):
    contact_profile = np.zeros(size, dtype=np.int32)

    chunk_rows, _ = contacts.chunks
    zeros = np.zeros(chunk_rows, dtype=np.int32)

    for chunk_beg in range(0, contacts.shape[0], chunk_rows):
        chunk_end = min(chunk_beg + chunk_rows, contacts.shape[0])
        ii, jj, cc = contacts[chunk_beg:chunk_end].astype(np.int32).T

        selector = (chain_ids[ii] == chain_ids[jj])
        ii = ii[selector]
        jj = jj[selector]
        cc = cc[selector]

        kk = np.abs(ii - jj)
        zz = zeros[:len(kk)]

        matrix = sparse.coo_matrix((cc, (kk, zz)), shape=(size, 1))
        contact_profile += matrix.toarray()[:, 0]

    return contact_profile


def fit_power_law(x, y):
    mask = (x > 0) & (y > 0)
    X = np.log(x[mask, None])
    y = np.log(y[mask])
    w = 1 / x[mask]
    linreg = LinearRegression()
    linreg.fit(X, y, w)
    return linreg.coef_[0], np.exp(linreg.intercept_)
