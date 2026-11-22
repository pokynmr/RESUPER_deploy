"""Tests for the ensemble writer: PDB round-trip and mmCIF output."""
from __future__ import annotations

import numpy as np
import pytest

from resuper.pdb import Residue, Structure, parse_pdb
from resuper.writer import superpose_ensemble, superpose_to_mean, write_cif, write_pdb


def _rot(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ])


def _base_coords():
    """Base xyz for every atom of the two residues (reference frame)."""
    return {
        ("A", 1): {"N": np.array([1.0, 2.0, 3.0]), "CA": np.array([2.0, 3.0, 4.0]),
                   "C": np.array([4.0, 1.0, 2.0]), "O": np.array([0.5, 0.5, 0.5]),
                   "CB": np.array([3.0, 3.0, 3.0])},
        ("A", 2): {"N": np.array([5.0, 1.0, 2.0]), "CA": np.array([6.0, 2.0, 3.0]),
                   "C": np.array([8.0, 0.0, 1.0]), "O": np.array([4.5, 4.5, 4.5]),
                   "CB": np.array([7.0, 7.0, 7.0])},
    }


def _make_structure() -> Structure:
    atoms = ("N", "CA", "C", "O", "CB")
    transforms = [(np.eye(3), np.zeros(3)),
                  (_rot(np.array([1.0, 2.0, 3.0]), 0.9), np.array([1.0, -2.0, 0.5])),
                  (_rot(np.array([-2.0, 1.0, 1.0]), 1.7), np.array([-0.5, 0.0, 1.5]))]
    residues = [Residue("A", 1, "ALA", 0), Residue("A", 2, "GLY", 1)]
    conformers = []
    for r, t in transforms:
        conf = {}
        for key, base in _base_coords().items():
            conf[key] = {a: base[a] @ r.T + t for a in base}
        conformers.append(conf)
    return Structure(residues=residues, conformers=conformers,
                     model_ids=[1, 2, 3])


def test_superpose_onto_first_model_restores_reference_frame():
    st = _make_structure()
    sup = superpose_ensemble(st, 0)
    ref = sup[0]
    for k in range(1, st.n_models):
        for key in ref:
            for a in ref[key]:
                np.testing.assert_allclose(sup[k][key][a], ref[key][a], atol=1e-6)


def test_write_pdb_round_trips_coordinates(tmp_path):
    st = _make_structure()
    sup = superpose_ensemble(st, 0)
    out = tmp_path / "ens.pdb"
    write_pdb(str(out), st, sup)
    reparsed = parse_pdb(str(out))
    assert reparsed.n_models == 3
    assert len(reparsed.residues) == 2
    for (chain, num), atoms_sup in sup[1].items():
        for a, xyz in atoms_sup.items():
            np.testing.assert_allclose(
                reparsed.conformers[1][(chain, num)][a], xyz, atol=1e-3)


def test_write_cif_emits_expected_atom_rows(tmp_path):
    st = _make_structure()
    sup = superpose_ensemble(st, 0)
    out = tmp_path / "ens.cif"
    write_cif(str(out), st, sup)
    text = out.read_text()
    assert "_atom_site.pdbx_PDB_model_num" in text
    atom_rows = [ln for ln in text.splitlines() if ln.startswith("ATOM ")]
    n_atoms = sum(len(c[(r.chain, r.num)]) for c in sup for r in st.residues)
    assert len(atom_rows) == n_atoms


def test_write_pdb_subset_selects_requested_models(tmp_path):
    st = _make_structure()
    sup = superpose_ensemble(st, 0)
    out = tmp_path / "subset.pdb"
    write_pdb(str(out), st, sup, model_indices=[1])     # second model only
    reparsed = parse_pdb(str(out))
    assert reparsed.n_models == 1
    for (chain, num), atoms_sup in sup[1].items():
        for a, xyz in atoms_sup.items():
            np.testing.assert_allclose(
                reparsed.conformers[0][(chain, num)][a], xyz, atol=1e-3)


def test_superpose_to_mean_aligns_ensemble():
    st = _make_structure()
    sup = superpose_to_mean(st)
    res = st.residues
    bb = ("N", "CA", "C")

    def _backbone(conf):
        return np.array([conf[(r.chain, r.num)][a] for r in res for a in bb])

    ref = _backbone(sup[0])
    for k in range(1, st.n_models):
        np.testing.assert_allclose(_backbone(sup[k]), ref, atol=1e-6)
