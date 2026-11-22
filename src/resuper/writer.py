"""Coordinate emission for RESUPER.

Write a multi-model ensemble as PDB or mmCIF, superposed onto a reference
conformer (by default the first model), optionally split by the detected
conformational states.

The superposition uses the backbone (N, CA, C) of every kept residue as the
fit set and rigidly applies the resulting transformation to all parsed atoms of
each conformer, so the written file contains every parsed atom in the
ensemble frame.  Only residues that survived :func:`resuper.pdb._complete`
(those with a full backbone in every conformer) are emitted; incomplete
residues are listed in ``Structure.dropped``.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from .geometry import _kabsch
from .pdb import BACKBONE, Structure, _atom_element


class _Transform:
    """A rigid transform (rotation R plus translation) derivable from Kabsch."""

    def __init__(self, r: np.ndarray, p_mean: np.ndarray, q_mean: np.ndarray):
        self.r = r
        self.p_mean = p_mean
        self.q_mean = q_mean

    def apply(self, xyz: np.ndarray) -> np.ndarray:
        return (xyz - self.p_mean) @ self.r.T + self.q_mean


def _fit_transform(conformer: Dict, reference: Dict,
                   residues) -> _Transform:
    """Rigid transform that superposes a conformer's backbone onto the
    reference's backbone."""
    p = np.array([conformer[(r.chain, r.num)][a] for r in residues for a in BACKBONE])
    q = np.array([reference[(r.chain, r.num)][a] for r in residues for a in BACKBONE])
    r = _kabsch(p, q)
    return _Transform(r, p.mean(axis=0), q.mean(axis=0))


def superpose_ensemble(structure: Structure, ref_index: int = 0
                       ) -> List[Dict]:
    """Return every conformer rigidly superposed onto the reference conformer.

    Each returned dict mirrors the input conformer layout (``(chain, resnum) ->
    {atom: xyz}``) with the transformed coordinates.  The reference conformer
    is unchanged (its transform is the identity up to its own fit).
    """
    ref = structure.conformers[ref_index]
    transforms = [_fit_transform(c, ref, structure.residues)
                  for c in structure.conformers]
    return [{key: {a: t.apply(xyz) for a, xyz in atoms.items()}
             for key, atoms in c.items()}
            for t, c in zip(transforms, structure.conformers)]


def superpose_to_mean(structure: Structure, iters: int = 4) -> List[Dict]:
    """Return every conformer superposed onto the ensemble-mean frame.

    Mirrors :func:`resuper.geometry.rmsd_to_mean`: iteratively Kabsch-fit each
    conformer (all its atoms) onto the mean backbone of the currently
    superposed set.  The result is the ensemble in the mean frame — a fixed
    reference rather than the first model.
    """
    residues = structure.residues
    cur = [{key: dict(atoms) for key, atoms in c.items()}
           for c in structure.conformers]
    for _ in range(iters):
        backs = [np.array([cw[(r.chain, r.num)][a]
                           for r in residues for a in BACKBONE]) for cw in cur]
        mean = np.mean(backs, axis=0)
        q_mean = mean.mean(axis=0)
        nxt = []
        for cc, b in zip(cur, backs):
            r = _kabsch(b, mean)
            p_mean = b.mean(axis=0)
            nxt.append({key: {a: (xyz - p_mean) @ r.T + q_mean
                              for a, xyz in atoms.items()}
                        for key, atoms in cc.items()})
        cur = nxt
    return cur


def _pdb_atom_line(serial: int, name: str, resname: str, chain: str, resnum: int,
                   xyz: np.ndarray, element: str) -> str:
    line = [" "] * 80
    line[0:6] = list("ATOM  ")
    line[6:11] = list(f"{serial:5d}")
    line[12:16] = list(f"{name:<4s}")
    line[17:20] = list(f"{resname:>3s}")
    line[21] = chain if chain else " "
    line[22:26] = list(f"{resnum:4d}")
    line[30:38] = list(f"{xyz[0]:8.3f}")
    line[38:46] = list(f"{xyz[1]:8.3f}")
    line[46:54] = list(f"{xyz[2]:8.3f}")
    line[54:60] = list("  1.00")
    line[60:66] = list("  0.00")
    line[76:78] = list(f"{element:>2s}")
    return "".join(line).rstrip()


def write_pdb(path: str, structure: Structure, superposed: List[Dict],
              model_indices: Optional[List[int]] = None) -> None:
    """Write a multi-MODEL PDB of the (superposed) conformers."""
    indices = list(range(len(structure.conformers))) if model_indices is None \
        else model_indices
    serial = 1
    with open(path, "w") as fh:
        for mi in indices:
            model_num = structure.model_ids[mi] \
                if mi < len(structure.model_ids) else mi + 1
            fh.write(f"MODEL     {model_num}\n")
            cm = superposed[mi]
            for r in structure.residues:
                atoms = cm.get((r.chain, r.num))
                if not atoms:
                    continue
                for name in sorted(atoms):
                    element = _atom_element(name) or name[:1]
                    fh.write(_pdb_atom_line(serial, name, r.name, r.chain, r.num,
                                            atoms[name], element) + "\n")
                    serial += 1
            fh.write("ENDMDL\n")
        fh.write("END\n")


_CIF_COLS = [
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_seq_id",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.pdbx_PDB_model_num",
]


def write_cif(path: str, structure: Structure, superposed: List[Dict],
              model_indices: Optional[List[int]] = None) -> None:
    """Write the (superposed) ensemble as an mmCIF file (multi-model)."""
    indices = list(range(len(structure.conformers))) if model_indices is None \
        else model_indices
    with open(path, "w") as fh:
        fh.write(f"data_{path.rsplit('/', 1)[-1].split('.')[0]}\n")
        fh.write("#\n# Submitted by RESUPER\n#\n")
        fh.write("loop_\n" + "".join(f"{c}\n" for c in _CIF_COLS))
        serial = 1
        for mi in indices:
            model_num = structure.model_ids[mi] \
                if mi < len(structure.model_ids) else mi + 1
            cm = superposed[mi]
            for r in structure.residues:
                atoms = cm.get((r.chain, r.num))
                if not atoms:
                    continue
                for name in sorted(atoms):
                    element = _atom_element(name) or name[:1]
                    x, y, z = atoms[name]
                    fh.write(
                        f"ATOM {serial} {element} {name} {r.name} {r.chain} "
                        f"{r.index + 1} {r.num} {r.name} {r.chain} {name} "
                        f"{x:.3f} {y:.3f} {z:.3f} 1.00 0.00 {model_num}\n")
                    serial += 1
        fh.write("#\n")
