"""Unit tests for the geometric primitives used by resuper."""
import numpy as np

from resuper.geometry import (dihedral, order_parameter, rmsd_sq,
                                rmsd_to_mean, rmsd_to_mean_k, superpose)


def _triangle():
    # Non-collinear 3-atom basis; Kabsch has full rotational DOF (no degeneracy).
    return np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


def _two_conf():
    # conf2 differs internally (one vertex moved) -> non-rigid conformational change.
    base = _triangle()
    c2 = base.copy()
    c2[1] = np.array([1.5, 0.0, 0.0])
    return np.stack([base, c2])


def _pairwise(ens):
    """RMSD between two conformers after Kabsch superposition onto the first."""
    sup = superpose(ens[1], ens[0])
    d2 = ((sup - ens[0]) ** 2).sum(axis=1).mean()
    return float(np.sqrt(d2))


def test_superpose_restores_rigid_translation():
    S = _triangle() + np.array([2.0, 5.0, 1.0])
    out = superpose(S, _triangle())
    assert np.allclose(out, _triangle(), atol=1e-9)


def test_rmsd_to_mean_half_pairwise():
    ens = _two_conf()
    p = _pairwise(ens)
    assert p > 0.01
    assert np.isclose(rmsd_to_mean(ens, ref_seq=ens[0]), 0.5 * p, atol=1e-9)


def test_rmsd_sq_equal_to_mean_for_n2():
    ens = _two_conf()
    # For n=2 the sqrt-mean-square and mean-of-RMSDs quantities coincide.
    assert np.isclose(rmsd_sq(ens, ref_seq=ens[0]),
                      rmsd_to_mean(ens, ref_seq=ens[0]), atol=1e-9)


def test_rmsd_to_mean_k_per_conformer():
    ens = _two_conf()
    p = _pairwise(ens)
    k = rmsd_to_mean_k(ens, ref=ens[0])
    assert k.shape == (2,)
    assert np.allclose(k, [0.5 * p, 0.5 * p], atol=1e-9)


def test_dihedral_known_geometry():
    p0 = np.array([1.0, 0, 0])
    p1 = np.array([0.0, 0, 0])
    p2 = np.array([0.0, 1.0, 0])
    p3 = np.array([0.0, 1.0, 1.0])
    d = dihedral(p0, p1, p2, p3)
    assert np.isclose(abs(d), np.pi / 2, atol=1e-9)


def test_order_parameter_ranges():
    assert np.isclose(order_parameter(np.zeros(10)), 1.0, atol=1e-9)
    # Uniformly spaced angles cancel to a near-zero order parameter.
    assert order_parameter(np.linspace(0, 2 * np.pi, 8, endpoint=False)) < 1e-9
