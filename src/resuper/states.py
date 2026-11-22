"""Multi-conformer state detection for RESUPER.

Given a bundle of conformers (multiple models of one structure), detect how
many distinct conformational states the ensemble adopts and report the
conformer-to-state assignment.

Descriptor: rotation-invariant triadistance PCA.  For every conformer the
upper-triangle of the CA-CA pairwise distance matrix is flattened into a
feature vector; the whole ensemble is mean-centered in that space and reduced
with PCA (via SVD).  Clustering is performed in the leading principal-component
space, which is rotation-invariant so the number of states is defined by shape
alone, not by arbitrary superposition.

Selection: the ensemble is a single state unless the silhouette of a clean
two-way split in the top-N principal-component space clears a calibrated floor;
otherwise the largest silhouette jump is used to pick the number of states.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage


@dataclass
class StateInfo:
    """A conformational state: its member conformers and how populated it is."""
    members: List[int]          # conformer numbers as the input file numbers them
    indices: List[int]          # 0-based positions in Structure.conformers
    population: float           # fraction of conformers in this state

    @property
    def size(self) -> int:
        return len(self.members)


@dataclass
class StatesReport:
    n_states: int
    labels: np.ndarray          # per-conformer state id (0-based)
    silhouette: Dict[int, float]
    states: List[StateInfo] = field(default_factory=list)


def conformer_numbers(structure) -> List[int]:
    """The number the input file gives each conformer, in conformer order.

    MODEL records win when they cover every conformer: a bundle numbered 21..40
    must not be reported as 1..20.  Bundles with no MODEL records (conformers
    split on a repeated atom) fall back to positional 1..N."""
    ids = list(structure.model_ids)
    if len(ids) == structure.n_models:
        return ids
    return list(range(1, structure.n_models + 1))


def triadistance_descriptor(structure) -> np.ndarray:
    """Rotation-invariant triadistance feature vector for every conformer.

    Returns a (N, n*(n-1)/2) float array: each row is the mean-centered
    upper-triangle CA-CA distance matrix of one conformer, as a flat vector.
    The descriptor is mean-centered but NOT z-scored (matching the scale the
    clustering was calibrated on).
    """
    res = structure.residues
    keys = [(r.chain, r.num) for r in res]
    n_ca = len(keys)
    iu = np.triu_indices(n_ca, k=1)
    rows = []
    for c in structure.conformers:
        ca = np.array([c[k]["CA"] for k in keys])
        d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
        rows.append(d[iu])
    mat = np.array(rows)
    return mat - mat.mean(axis=0, keepdims=True)


def _principal_component_space(descriptor: np.ndarray, ncomp: int) -> np.ndarray:
    """Scores of the leading principal components (rotation-invariant space)."""
    if ncomp >= descriptor.shape[1]:
        ncomp = descriptor.shape[1] - 1
    ncomp = max(1, ncomp)
    u, s, _ = np.linalg.svd(descriptor, full_matrices=False)
    # Scores = U * S, first ncomp columns (the variance-bearing axes).
    return (u * s)[:, :ncomp]


def _silhouette(xc: np.ndarray, lab: np.ndarray) -> float:
    """Mean silhouette over all conformers in the given feature space."""
    labs = np.unique(lab)
    if labs.size <= 1:
        return 0.0
    dmat = np.linalg.norm(xc[:, None] - xc[None], axis=2)
    s = np.empty(xc.shape[0])
    for i in range(xc.shape[0]):
        same = lab == lab[i]
        a = dmat[i, same].mean() if same.sum() > 1 else 0.0
        neighbors = [
            dmat[i, lab == lb].mean()
            for lb in labs if lb != lab[i] and (lab == lb).any()
        ]
        b = min(neighbors) if neighbors else 0.0
        denom = max(a, b)
        s[i] = (b - a) / denom if denom > 0 else 0.0
    return float(s.mean())


def _cluster_space(xc: np.ndarray, k: int, min_state: int) -> Optional[np.ndarray]:
    """Ward clustering to `k` groups; return labels or None if no group is too
    small (a split that isolates fewer than `min_state` conformers is not a
    valid two-state decision)."""
    if xc.shape[0] < 2 * min_state or k > xc.shape[0] // min_state:
        return None
    z = linkage(xc, method="ward")
    lab = fcluster(z, k, criterion="maxclust").astype(int)
    # fcluster yields 1-indexed labels; check every used group's size.
    sizes = np.bincount(lab)[np.unique(lab)]
    if sizes.min() < min_state:
        return None
    return lab


def detect_states(structure, ncomp: int = 3, floor: float = 0.60,
                  min_state: int = 3) -> StatesReport:
    """Detect the number of conformational states in a conformer bundle.

    Returns a :class:`StatesReport` with the number of states, the
    conformer-to-state assignment, silhouette scores per split, and per-state
    population breakdown.

    Parameters
    ----------
    structure : Structure
        Parsed multi-model bundle.
    ncomp : int
        Number of leading principal components retained (default 3).
    floor : float
        Minimum two-state silhouette for the ensemble to be called two states.
    min_state : int
        Smallest group size for a split to count as a valid two-state decision.
    """
    n = structure.n_models
    if n < 2:
        # Nothing to split: trivially a single state.
        return StatesReport(
            n_states=1,
            labels=np.zeros(n, dtype=int),
            silhouette={},
            states=[StateInfo(conformer_numbers(structure), list(range(n)), 1.0)],
        )

    xc = _principal_component_space(triadistance_descriptor(structure), ncomp)

    # Silhouette + labels for every valid split size, from 2 up to the largest
    # where each group still holds at least `min_state` conformers.
    silhouette: Dict[int, float] = {1: 0.0}
    _labels: Dict[int, np.ndarray] = {1: np.zeros(n, dtype=int)}
    for k in range(2, n // min_state + 1):
        lab = _cluster_space(xc, k, min_state)
        if lab is None:
            continue
        silhouette[k] = _silhouette(xc, lab)
        _labels[k] = lab

    # Choose the number of states.
    best_k = 1
    best_jump = -np.inf
    for k in sorted(silhouette):
        if k == 1:
            continue
        jump = silhouette[k] - silhouette.get(k - 1, 0.0)
        if jump > best_jump:
            best_jump, best_k = jump, k

    # Two states only if the clean two-way split clears the silhouette floor.
    n_states = best_k if (best_k >= 2 and silhouette[best_k] >= floor) else 1

    # Canonicalize labels for the *chosen* number of states so the split
    # ordering is stable across calls (the all-in-one case has no split labels).
    chosen_labels = _labels.get(n_states, np.zeros(n, dtype=int))
    ordered = np.unique(chosen_labels)
    canon = np.searchsorted(ordered, chosen_labels)

    return StatesReport(
        n_states=n_states,
        labels=canon,
        silhouette=silhouette,
        states=_state_infos(canon, conformer_numbers(structure)),
    )


def _state_infos(labels: np.ndarray, numbers: List[int]) -> List[StateInfo]:
    """Group conformers by their (canonical, 0-based) state id."""
    infos: List[StateInfo] = []
    n = len(labels)
    for state in range(int(labels.max()) + 1):
        idx = [i for i in range(n) if labels[i] == state]
        infos.append(StateInfo([numbers[i] for i in idx], idx, len(idx) / n))
    return infos
