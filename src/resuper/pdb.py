"""Multi-model PDB parsing for RESUPER.

RESUPER operates on a bundle of conformers (multiple MODEL records) of the same
sequence. This module parses such a PDB into ordered residues and per-conformer
atomic coordinates, and provides torsion-angle helpers.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

# Backbone atoms used for superposition/RMSD.
BACKBONE = ("N", "CA", "C")
# Carbonyl "C" in PDB is the C' (prime) backbone atom.

# The 20 standard amino-acid names. A residue that is a HETATM record OR carries
# a non-standard (not in this set) residue name is not part of the protein chain
# and is skipped when building the residue sequence, which shows up as an
# internal gap in the reported optimal range (e.g. 1A3P's ABA at res6/20).
_STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@dataclass
class Residue:
    chain: str
    num: int            # residue number as int (insertion codes ignored)
    name: str
    index: int          # 0-based position in the sequence of residues


@dataclass
class Structure:
    residues: List[Residue]                    # sequence order
    conformers: List[Dict[str, np.ndarray]]    # per model: (chain,resnum) -> {atom: xyz}
    model_ids: List[int] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)   # report lines (format_report)
    dropped: List[Tuple[str, int]] = field(default_factory=list)  # incomplete residues
    few_core: bool = False      # fewer core atoms than the minimal domain size

    @property
    def n_models(self) -> int:
        return len(self.conformers)


def _atom_element(name: str) -> str:
    # PDB columns 13-14 give element; fall back to name heuristics.
    return name[:1] if name else ""


def parse_pdb(path: str) -> Structure:
    residues: List[Residue] = []
    conformers: List[Dict[str, np.ndarray]] = []
    model_ids: List[int] = []

    cur: Optional[Dict[str, np.ndarray]] = None
    seen: Dict[Tuple[str, int], Residue] = {}
    order: List[Tuple[str, int]] = []

    def ensure_residue(chain, resnum, resname):
        key = (chain, resnum)
        if key not in seen:
            res = Residue(chain=chain, num=resnum, name=resname.upper(), index=len(order))
            seen[key] = res
            residues.append(res)
            order.append(key)
        return key

    with open(path) as fh:
        for line in fh:
            rec = line[0:6].strip().upper()
            if rec == "MODEL":
                tokens = line.split()
                model_ids.append(int(tokens[1]) if len(tokens) > 1 else len(model_ids) + 1)
                cur = {}
                continue
            if rec == "ENDMDL":
                if cur is not None:
                    conformers.append(cur)
                    cur = None
                continue
            at_het = line.startswith("HETATM")
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if cur is None:
                cur = {}
            name = line[12:16].strip()
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                resnum = int(line[22:26].strip() or 0)
            except ValueError:
                continue                # unreadable record: skip it
            if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)):
                continue                # NaN/inf coordinate: treat the atom as missing
            chain = line[21].strip() or "A"
            resname = line[17:20].strip()
            # HETATM and non-standard residues are not part of the protein
            # chain (seen as an interior gap in the reported range, e.g. 1A3P
            # ABA).
            if at_het or resname not in _STANDARD_RESIDUES:
                continue
            key = ensure_residue(chain, resnum, resname)
            atoms = cur.setdefault(key, {})
            if name in atoms:
                if not model_ids and line[16] == " ":
                    # Without MODEL records a repeated atom starts the next model.
                    conformers.append(cur)
                    cur = {key: {}}
                    atoms = cur[key]
                else:
                    continue            # alternate location or duplicate: keep the first
            atoms[name] = np.array([x, y, z], dtype=float)

    if cur is not None:
        conformers.append(cur)
    return _complete(Structure(residues=residues, conformers=conformers, model_ids=model_ids))


def _complete(st: Structure) -> Structure:
    """Drop residues that lack N, CA or C in any conformer (they cannot be superimposed
    or measured) and re-index the rest; their (chain, number) go to `st.dropped`."""
    keep = [r for r in st.residues
            if all(all(a in c.get((r.chain, r.num), {}) for a in BACKBONE)
                   for c in st.conformers)]
    if len(keep) == len(st.residues):
        return st
    st.dropped = [(r.chain, r.num) for r in st.residues if r not in keep]
    for i, r in enumerate(keep):
        r.index = i
    st.residues = keep
    return st


def subset_conformers(st: Structure, indices) -> Structure:
    """A Structure holding only the conformers at `indices`, sequence unchanged.

    Used to analyse one conformational state as an ensemble of its own. Nothing
    needs re-superposing: Stage A works on inter-conformer distances and Stage B
    fits its own transform per range, so both are independent of the frame the
    file happens to sit in. Nor does anything need re-completing -- a residue
    that carries N, CA and C in every conformer of the bundle certainly does so
    in a subset of them.

    `warnings` and `few_core` are the full bundle's Stage A verdict, so they do
    not carry over: the subset earns its own.
    """
    indices = list(indices)
    model_ids = ([st.model_ids[i] for i in indices]
                 if len(st.model_ids) == st.n_models else [])
    return Structure(residues=st.residues,
                     conformers=[st.conformers[i] for i in indices],
                     model_ids=model_ids, warnings=[], dropped=list(st.dropped),
                     few_core=False)

