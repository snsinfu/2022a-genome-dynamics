"""
This module defines classes to access simulation trajectory file.
"""

import json

from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from types import SimpleNamespace

import h5py
import numpy as np


_GROUP_PHASE_PARENT = "snapshots"
_GROUP_METADATA = "metadata"

_DATASET_CONFIG = "config"
_DATASET_PARTICLE_TYPES = "particle_types"
_DATASET_AB_FACTORS = "ab_factors"
_DATASET_CHROMOSOME_RANGES = "chromosome_ranges"
_DATASET_CENTROMERE_RANGES = "centromere_ranges"
_DATASET_NUCLEOLUS_RANGES = "nucleolus_ranges"
_DATASET_STEPS = ".steps"
_DATASET_CONTEXT = "context"
_DATASET_POSITIONS = "positions"
_DATASET_CONTACT_MAP = "contact_map"


@dataclass
class Span:
    """
    Named interval in a flattened dataset.
    """
    name: str
    start: int
    end: int


@dataclass
class Chromosome(Span):
    """
    Chromosome interval in a flattened dataset.
    """
    centromere_start: int
    centromere_end: int


@dataclass
class Metadata:
    """
    Simuation metadata.
    """
    config: SimpleNamespace
    particle_types: np.ndarray
    particle_types_enum: Dict[str, int]
    ab_factors: np.ndarray
    chromosomes: List[Chromosome]
    nucleoli: List[Span]


class Snapshot:
    """
    Single snapshot sample.
    """
    def __init__(self, node: h5py.Group, step: int):
        self._node = node
        self._step = step

    @property
    def step(self) -> int:
        """
        Returns the simulation step number from which the snapshot is taken.
        """
        return self._step

    @property
    def positions(self) -> np.ndarray:
        """
        Returns an array of coordinate values of the particles.
        """
        return self._node[_DATASET_POSITIONS][:]

    @property
    def context(self) -> SimpleNamespace:
        """
        Returns simulation context and statistics.
        """
        return SimpleNamespace(**json.loads(self._node[_DATASET_CONTEXT][()]))

    @property
    def contact_map_dataset(self) -> Optional[h5py.Dataset]:
        """
        Returns a dataset containing (i,j,v)-formatted contact map, if any.
        """
        if _DATASET_CONTACT_MAP in self._node:
            return self._node[_DATASET_CONTACT_MAP]
        return None


class Phase:
    """
    Simulation phase.
    """
    def __init__(self, branch: h5py.Group):
        self._branch = branch
        self._metadata = _make_phase_metadata(self._branch)
        self._snapshots = _make_snapshots(self._branch)

    @property
    def metadata(self) -> Metadata:
        """
        Returns the `Metadata` used in the phase.
        """
        return self._metadata

    @property
    def snapshots(self) -> List[Snapshot]:
        """
        Returns a list of `Snapshot` samples in the phase.
        """
        return self._snapshots


def _make_snapshots(branch):
    snapshots = []
    for step in branch[_DATASET_STEPS]:
        snapshots.append(Snapshot(branch[step], int(step)))
    return snapshots


def _make_phase_metadata(branch):
    default_metadata = branch.file[_GROUP_METADATA]
    if _GROUP_METADATA in branch:
        metadata = branch[_GROUP_METADATA]
    else:
        metadata = default_metadata

    def dataset(key):
        if key in metadata:
            return metadata[key]
        return default_metadata[key]

    # Global configurations.
    config = SimpleNamespace(**json.loads(dataset(_DATASET_CONFIG)[()]))

    # Particle parameters.
    particle_types = dataset(_DATASET_PARTICLE_TYPES)
    particle_types_enum = dict(h5py.check_dtype(enum=particle_types.dtype))
    particle_types = particle_types[:]

    ab_factors = dataset(_DATASET_AB_FACTORS)[:]

    # Chromosomes.
    chromosome_ranges = dataset(_DATASET_CHROMOSOME_RANGES)
    chromosome_keys = json.loads(chromosome_ranges.attrs["keys"])
    chromosome_ranges = chromosome_ranges[:]

    centromere_ranges = dataset(_DATASET_CENTROMERE_RANGES)
    centromere_keys = json.loads(centromere_ranges.attrs["keys"])
    centromere_ranges = centromere_ranges[:]

    chromosomes = []
    for name in chromosome_keys:
        start, end = chromosome_ranges[chromosome_keys[name]]
        cen_start, cen_end = centromere_ranges[centromere_keys[name]]
        chrom = Chromosome(name, start, end, cen_start, cen_end)
        chromosomes.append(chrom)

    # Nucleolus.
    nucleolus_ranges = dataset(_DATASET_NUCLEOLUS_RANGES)
    nucleolus_keys = json.loads(nucleolus_ranges.attrs["keys"])
    nucleolus_ranges = nucleolus_ranges[:]

    nucleoli = []
    for name in nucleolus_keys:
        start, end = nucleolus_ranges[nucleolus_keys[name]]
        span = Span(name, start, end)
        nucleoli.append(span)

    return Metadata(
        config, particle_types, particle_types_enum, ab_factors, chromosomes, nucleoli
    )


class Store:
    """
    Simulation trajectory file.
    """
    def __init__(self, path: str, mode: str = "r"):
        self._file = h5py.File(path, mode)

    def __enter__(self):
        self._file.__enter__()
        return self

    def __exit__(self, *args):
        self._file.__exit__(*args)

    def phase(self, name: str) -> Phase:
        """
        Returns the simulation phase.
        """
        return Phase(self._file[_GROUP_PHASE_PARENT][name])
