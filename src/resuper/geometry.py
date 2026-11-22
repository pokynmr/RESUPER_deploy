"""Geometric primitives for RESUPER: dihedral angles, Kabsch superposition,
mean-coordinate RMSD, and angular order parameters."""
from __future__ import annotations

import numpy as np


def dihedral(p0, p1, p2, p3) -> float:
    """Signed dihedral angle (radians) around the p1-p2 bond, in [-pi, pi]."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - b1n * np.dot(b0, b1n)
    w = b2 - b1n * np.dot(b2, b1n)
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return np.arctan2(y, x)


def _kabsch(P, Q):
    """Rotation matrix R (3x3) superposing P onto Q (both (M,3), same order)."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Dm = np.diag([1.0, 1.0, d])
    R = Vt.T @ Dm @ U.T
    return R


def superpose(S, ref):
    """Superpose array S (M,3) onto ref (M,3); returns rotated S."""
    R = _kabsch(S, ref)
    return (S - S.mean(axis=0)) @ R.T + ref.mean(axis=0)


def rmsd_to_mean(ens: np.ndarray, ref_seq=None, iters: int = 1) -> float:
    """Mean-coordinate RMSD of an ensemble, as RESUPER reports it.

    `ens` shape (N, M, 3) = (conformers, atoms, xyz). Superpose all conformers
    onto the reference (`ref_seq`, default the mean) by Kabsch, compute the mean
    coordinates of the superposed ensemble, and return the *average over
    conformers* of each conformer's RMSD to that mean.
    """
    N, M, _ = ens.shape
    mean = ens.mean(axis=0)
    ref = ref_seq if ref_seq is not None else mean
    for _ in range(iters):
        sup = np.array([superpose(ens[k], ref) for k in range(N)])
        mean = sup.mean(axis=0)
        if ref_seq is None:
            ref = mean
    d2 = ((sup - mean) ** 2).sum(axis=2).mean(axis=1)
    return float(np.sqrt(d2).mean())


def rmsd_to_mean_converged(ens: np.ndarray, tol: float = 1e-10, max_iter: int = 200) -> float:
    """Average per-conformer RMSD to the mean after iterative superposition.

    Conformers are first superposed onto the first conformer; the mean is then
    used as the new reference and the superposition repeated until the mean
    structure stops changing (generalized Procrustes fit).
    """
    N = ens.shape[0]
    sup = np.array([superpose(ens[k], ens[0]) for k in range(N)])
    for _ in range(max_iter):
        mean = sup.mean(axis=0)
        new = np.array([superpose(ens[k], mean) for k in range(N)])
        done = np.abs(new - sup).max() < tol
        sup = new
        if done:
            break
    mean = sup.mean(axis=0)
    d2 = ((sup - mean) ** 2).sum(axis=2).mean(axis=1)
    return float(np.sqrt(d2).mean())


def rmsd_sq(ens: np.ndarray, ref_seq=None, iters: int = 1) -> float:
    """Sqrt-of-mean-square RMSD to the mean (the internal `r` of refinement step 1).

    Identical superposition/mean construction as `rmsd_to_mean`, but the returned
    value is the global sqrt of mean squared distances rather than the average of
    per-conformer RMSDs. The reference mode uses this quantity for the
    removal-decision thresholds while reporting `rmsd_to_mean`.
    """
    N, M, _ = ens.shape
    mean = ens.mean(axis=0)
    ref = ref_seq if ref_seq is not None else mean
    for _ in range(iters):
        sup = np.array([superpose(ens[k], ref) for k in range(N)])
        mean = sup.mean(axis=0)
        if ref_seq is None:
            ref = mean
    d2 = ((sup - mean) ** 2).sum(axis=2)
    return float(np.sqrt(d2.mean()))


def rmsd_to_mean_k(ens: np.ndarray, ref=None) -> np.ndarray:
    """Per-conformer RMSD to the mean (N,) after superposition onto `ref`."""
    N, M, _ = ens.shape
    ref = ref if ref is not None else ens.mean(axis=0)
    sup = np.array([superpose(ens[k], ref) for k in range(N)])
    mean = sup.mean(axis=0)
    d2 = ((sup - mean) ** 2).sum(axis=2).mean(axis=1)
    return np.sqrt(d2)


def order_parameter(angles: np.ndarray) -> float:
    """Angular order parameter S = |mean_k exp(i theta_k)| in [0,1]."""
    z = np.exp(1j * np.asarray(angles, dtype=float))
    return float(abs(z.mean()))
