"""RESUPER core: ordered-region identification and superposition ranges.

RESUPER builds on the two-stage framework of Kirchner & Guentert, BMC
Bioinformatics 2011 (12:170): (A) domain identification via angular order
parameters and distance-variance clustering; (B) per-domain residue-range
refinement via an iterative gap-weighted RMSD-decrease procedure. The
`default` method replaces the opaque/tuned parts of the published framework
with decisions that have a geometric meaning, and is the method this package
provides. The `reference` method reproduces the CYRANGE 2.0 report and is kept
only for side-by-side comparison.

Constants are centralised below and in Params.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, to_tree
from scipy.spatial.distance import squareform

from .geometry import dihedral, rmsd_to_mean, rmsd_to_mean_converged, rmsd_sq, superpose, order_parameter
from .pdb import BACKBONE, Structure, subset_conformers

# Side-chain torsions per residue type (IUPAC chi1..chi5). Backbone phi/psi
# are added for every residue; "-"/"+" prefixes name the previous/next residue.
_CHI: Dict[str, Tuple[Tuple[str, ...], ...]] = {
    "ARG": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"),
            ("CG", "CD", "NE", "CZ"), ("CD", "NE", "CZ", "NH1")),
    "ASN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")),
    "ASP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")),
    "CYS": (("N", "CA", "CB", "SG"),),
    "GLN": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")),
    "GLU": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")),
    "HIS": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")),
    "ILE": (("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")),
    "LEU": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "LYS": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"),
            ("CG", "CD", "CE", "NZ")),
    "MET": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")),
    "PHE": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "PRO": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")),
    "SER": (("N", "CA", "CB", "OG"),),
    "THR": (("N", "CA", "CB", "OG1"),),
    "TRP": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "TYR": (("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")),
    "VAL": (("N", "CA", "CB", "CG1"),),
}


# Rotatable hydrogen groups (methyl, hydroxyl, thiol, amine): one dihedral each, defined
# by one of the hydrogens. Used only by the default method when DEFAULT_TORSIONS == "chi_h".
_H_ROTORS: Dict[str, Tuple[Tuple[str, ...], ...]] = {
    "ALA": (("N", "CA", "CB", "HB1"),),
    "VAL": (("CA", "CB", "CG1", "HG11"), ("CA", "CB", "CG2", "HG21")),
    "LEU": (("CB", "CG", "CD1", "HD11"), ("CB", "CG", "CD2", "HD21")),
    "ILE": (("CA", "CB", "CG2", "HG21"), ("CB", "CG1", "CD1", "HD11")),
    "THR": (("CA", "CB", "OG1", "HG1"), ("CA", "CB", "CG2", "HG21")),
    "SER": (("CA", "CB", "OG", "HG"),),
    "CYS": (("CA", "CB", "SG", "HG"),),
    "MET": (("CG", "SD", "CE", "HE1"),),
    "TYR": (("CE1", "CZ", "OH", "HH"),),
    "LYS": (("CD", "CE", "NZ", "HZ1"),),
}


@dataclass
class Params:
    absdecrease: float = 1.6     # delta_abs (Angstrom)
    buffer: int = 3              # m: domain-boundary extension (residues)
    decrease: float = 1.2        # delta (relative-decrease base)
    gap: float = 0.4             # gamma: gap-opening penalty
    rmsd_threshold: float = 0.1  # u (usage inferred; reserved)
    gap_width: int = 3           # g: small gaps (< g residues) are refilled
    lib: Optional[str] = None
    min_cluster: int = 8         # mu: minimal cluster/domain size
    residue_range: Optional[str] = None
    # "default": RESUPER's own method (published Q_i cutoff, explained-RMSD
    # domain splits, consistent domain ranges, uniform N/CA/C atom set in Stage B).
    # "reference": reproduce the CYRANGE 2.0 report for side-by-side comparison.
    method: str = "default"


@dataclass
class Domain:
    ranges: List[Tuple[int, int]]          # inclusive 1-based residue ranges
    rmsd: float
    n_residues: int
    n_gaps: int
    chains: List[str] = field(default_factory=list)  # per-range chain prefix


def _atom(structure, conf, chain, resnum, name):
    return conf[(chain, resnum)][name]


def _torsion_defs(name: str, h_rotors: bool = False):
    defs = [("phi", ("-C", "N", "CA", "C")), ("psi", ("N", "CA", "C", "+N"))]
    for k, seq in enumerate(_CHI.get(name, ())):
        defs.append((f"chi{k + 1}", seq))
    if h_rotors:
        for k, seq in enumerate(_H_ROTORS.get(name, ())):
            defs.append((f"h{k + 1}", seq))
    return defs


def residue_torsion_angles(structure: Structure, h_rotors: bool = False):
    """Per residue: dict torsion-name -> array of angle values across conformers.

    A torsion is reported only if all four atoms are present in every conformer.
    Angles are in radians.
    """
    res = structure.residues
    values: List[Dict[str, np.ndarray]] = [dict() for _ in res]
    for idx in range(len(res)):
        r = res[idx]
        for name, seq in _torsion_defs(r.name, h_rotors):
            pts = []
            ok = True
            for c in structure.conformers:
                coords = []
                for at in seq:
                    if at[0] in "-+":
                        j = idx - 1 if at[0] == "-" else idx + 1
                        if not (0 <= j < len(res)) or res[j].chain != r.chain:
                            ok = False
                            break
                        key = (res[j].chain, res[j].num)
                        at = at[1:]
                    else:
                        key = (r.chain, r.num)
                    a = c.get(key, {}).get(at)
                    if a is None:
                        ok = False
                        break
                    coords.append(a)
                if not ok:
                    break
                pts.append(coords)
            if not ok or len(pts) != structure.n_models:
                continue
            arr = np.array(pts)  # (N,4,3)
            theta = np.array([dihedral(*arr[k]) for k in range(arr.shape[0])])
            values[idx][name] = theta
    return values


def order_parameters(structure: Structure, h_rotors: bool = False) -> List[Dict[str, float]]:
    angles = residue_torsion_angles(structure, h_rotors)
    return [{k: order_parameter(v) for k, v in d.items()} for d in angles]


def _pick_cutoff(torsion_s: np.ndarray) -> float:
    """S_cut chosen to maximize Q_i (rank/range-normalized objective, paper Fig).

    Implemented from the prose: maximize over torsion i of
      Q_i = ((S_i - S_min)/(S_max - S_min)) * (r_i / s)
    where r_i is the rank of torsion i when S values are sorted ascending and s is
    the total number of torsions. (Exact form is an image in the paper; tuned
    against the reference implementation.)
    """
    if torsion_s.size == 0:
        return 0.0
    s = torsion_s.size
    smin, smax = torsion_s.min(), torsion_s.max()
    if smax - smin < 1e-12:
        return smin
    order = np.argsort(torsion_s, kind="stable")   # ascending; rank = position
    ranks = np.empty_like(order)
    ranks[order] = np.arange(s) + 1                # 1..s
    q = ((torsion_s - smin) / (smax - smin)) * (ranks / s)
    return float(torsion_s[int(np.argmax(q))])


def core_residue_flags(structure: Structure, s_cut: float) -> np.ndarray:
    """Boolean per residue: CA is a core atom (any torsion with S > s_cut)."""
    o = order_parameters(structure)
    flags = np.zeros(len(structure.residues), dtype=bool)
    for idx, d in enumerate(o):
        if any(v > s_cut for v in d.values()):
            flags[idx] = True
    return flags


def _c_terminal_mask(structure: Structure) -> np.ndarray:
    """Boolean per residue: the C-terminal residue as the reference sees it.

    Only the last residue of the whole structure is C-terminal. Earlier chain
    ends (with or without OXT) are ordinary residues: the reference reads all
    chains as one sequence (1E91 chain A end A85 keeps N/CA/C; 2I94 A189 too).
    """
    n = len(structure.residues)
    mask = np.zeros(n, dtype=bool)
    if n:
        mask[-1] = True
    return mask


def _backbone_matrix(structure, conformers, rng: Sequence[int],
                     exclude_cterm: bool = True):
    """Return ens (N,m,3) backbone coordinates for residues in rng, plus atom labels.

    `rng` is a list of 0-based residue indices. Each residue contributes its
    N, CA, C atoms, except a chain C-terminal residue, which contributes only
    its N (measured with frozen-range reference runs: the terminal's CA
    and C never enter the superposition/RMSD atom set). `exclude_cterm=False`
    keeps the full N/CA/C of the C-terminal.
    """
    N = len(conformers)
    cterm = _c_terminal_mask(structure)
    mats = []
    for idx in rng:
        names = ("N",) if (exclude_cterm and cterm[idx]) else BACKBONE
        for aname in names:
            mats.append((idx, aname))
    ens = np.zeros((N, len(mats), 3))
    for k, c in enumerate(conformers):
        for j, (idx, aname) in enumerate(mats):
            r = structure.residues[idx]
            ens[k, j] = c[(r.chain, r.num)][aname]
    return ens, mats


def _gap_count(sel: np.ndarray) -> int:
    """Number of gaps (maximal internal runs of non-selected residues)."""
    n = sel.size
    gaps = 0
    i = 0
    while i < n:
        if not sel[i]:
            j = i
            while j < n and not sel[j]:
                j += 1
            # internal gap: selected residue before i and after j-1
            if 0 < i and j < n:
                gaps += 1
            i = j
        else:
            i += 1
    return gaps


def _segments(sel: np.ndarray) -> List[Tuple[int, int]]:
    segs = []
    i = 0
    n = sel.size
    while i < n:
        if sel[i]:
            j = i
            while j < n and sel[j]:
                j += 1
            segs.append((i, j - 1))
            i = j
        else:
            i += 1
    return segs


def _number_segments(sel: np.ndarray, residues) -> Tuple[List[Tuple[int, int]], int]:
    """Split a selected mask into PDB-residue-number segments.

    Reports contiguous runs of selected residues split at any residue-number
    discontinuity. Non-standard residues are excluded from the chain at parse
    time (e.g. 1A3P's ABA at res20), so a jump in residue numbers corresponds to
    an excluded residue and is seen by the reference implementation as an internal gap (e.g. the
    the reference reports 1A3P as `7..19, 21..47`, not a contiguous `7..48`). The gap
    count is the number of internal number runs not selected.
    """
    segs = _chain_segments(sel, residues)
    # A gap is a maximal internally-missing number run, i.e. a boundary between
    # two reported segments. (Leading/trailing exclusions are not gaps.)
    return [(lo, hi) for _, lo, hi in segs], len(segs) - 1 if segs else 0


def _chain_segments(sel: np.ndarray, residues) -> List[Tuple[str, int, int]]:
    """Selected runs as (chain, lo, hi), split at number jumps and chain changes."""
    n = len(residues)
    segs: List[Tuple[str, int, int]] = []
    i = 0
    while i < n:
        if not sel[i]:
            i += 1
            continue
        chain = residues[i].chain
        lo = residues[i].num
        prev = lo
        i += 1
        while i < n and sel[i]:
            num = residues[i].num
            if num != prev + 1 or residues[i].chain != chain:
                segs.append((chain, lo, prev))
                chain = residues[i].chain
                lo = num
            prev = num
            i += 1
        segs.append((chain, lo, prev))
    return segs


class RangeRefiner:
    """Stage B: iterative residue-range refinement (7 steps, paper section 4)."""

    def __init__(self, structure: Structure, params: Params, residue_indices: Sequence[int]):
        self.structure = structure
        self.params = params
        self.idx = list(residue_indices)         # domain residues (0-based)
        self.conformers = structure.conformers
        self._default = params.method == "default"
        self._and_gate = self._default and DEFAULT_AND_GATE
        # The reference mode uses only the N atom of the C-terminal residue; the
        # default method treats every residue alike (N, CA, C).
        self._cterm = (np.zeros(len(structure.residues), dtype=bool) if self._default
                       else _c_terminal_mask(structure))

    def _is_cterm(self, pos: int) -> bool:
        """Whether the selected residue at domain-relative `pos` is a chain C-terminal."""
        return bool(self._cterm[self.idx[pos]])

    def _rmsd(self, sel: np.ndarray) -> Tuple[float, List[int]]:
        rng = [i for i in range(self.idx[0], self.idx[-1] + 1)
               if sel[i - self.idx[0]]]
        ens, _ = _backbone_matrix(self.structure, self.conformers, rng,
                                  exclude_cterm=not self._default)
        ref = ens[0]
        # The RMSD value r of step 1 is the quantity that is also reported, i.e.
        # the average over conformers of the per-conformer RMSD to the mean.
        r = rmsd_to_mean(ens, ref_seq=ref, iters=1)
        return r, rng

    def _reported_rmsd(self, sel: np.ndarray) -> float:
        rng = [i for i in range(self.idx[0], self.idx[-1] + 1)
               if sel[i - self.idx[0]]]
        ens, _ = _backbone_matrix(self.structure, self.conformers, rng,
                                  exclude_cterm=not self._default)
        return rmsd_to_mean_converged(ens)

    def _displacements(self, sel: np.ndarray) -> Dict[int, float]:
        rng = [i for i in range(self.idx[0], self.idx[-1] + 1)
               if sel[i - self.idx[0]]]
        ens, mats = _backbone_matrix(self.structure, self.conformers, rng,
                                  exclude_cterm=not self._default)
        ref = ens[0]
        sup = np.array([superpose(ens[k], ref) for k in range(ens.shape[0])])
        mean = sup.mean(axis=0)
        # Average displacement of a residue = mean over conformers and its
        # backbone atoms of the distance to the mean position. Residues excluded
        # from the superposition atom set (C-terminal) have no displacement (0).
        disp: Dict[int, float] = {idx: 0.0 for idx in rng}
        cols: Dict[int, List[int]] = {}
        for j, (idx, aname) in enumerate(mats):
            cols.setdefault(idx, []).append(j)
        for idx, js in cols.items():
            dist = np.linalg.norm(sup[:, js] - mean[js], axis=2)
            disp[idx] = float(dist.mean())
        return disp

    def _opens_gap(self, sel: np.ndarray, pos: int) -> bool:
        """Removing the selected residue at `pos` opens a new gap."""
        cur = sel.copy()
        cur[pos] = False
        return _gap_count(cur) > _gap_count(sel)

    def _rmsd_remove(self, sel: np.ndarray, pos: int) -> float:
        cur = sel.copy()
        cur[pos] = False
        r_new, _ = self._rmsd(cur)
        return r_new

    def _decrease(self, sel: np.ndarray, pos: int, r: float) -> float:
        r_new = self._rmsd_remove(sel, pos)
        dr = r - r_new
        if self._opens_gap(sel, pos):
            dr = self.params.gap * dr
        return dr

    def _condition_remove(self, sel: np.ndarray, pos: int, r: float, dr: float) -> bool:
        """Step 5 condition: dr >= delta_abs*n/N and dr/r >= delta_rel*n/N.

        n/N = atoms of the removed residue / atoms of all selected residues
        (each residue contributes len(BACKBONE) atoms), so n/N = 1/M.
        """
        p = self.params
        M = int(sel.sum())
        # n/N = atoms of the removed residue / atoms of all selected residues
        # (a chain C-terminal contributes only its N atom, see _backbone_matrix).
        n_atoms = lambda q: 1 if self._is_cterm(q) else len(BACKBONE)
        N_atoms = sum(n_atoms(q) for q in range(len(self.idx)) if sel[q])
        n_over_N = n_atoms(pos) / N_atoms
        delta_rel = p.decrease + 3.0 / M   # delta^rel = delta + 3.0/M
        # The reference removes a residue when its gap-weighted RMSD decrease is
        # large enough in absolute terms (>= delta_abs*n/N) OR in relative terms
        # (dr/r >= delta_rel*n/N) -- the two gates are alternative, not
        # conjunctive (verified against the reference implementation: 1A67's second gap opens via
        # the absolute gate while the relative gate just misses; an AND of the
        # two over-tightens and fails the boundary cases).
        abs_ok = dr >= p.absdecrease * n_over_N
        rel_ok = dr / r >= delta_rel * n_over_N
        # DEFAULT_AND_GATE switches the default method to the paper text's AND.
        if dr >= 0.0 and ((abs_ok and rel_ok) if self._and_gate else (abs_ok or rel_ok)):
            sel[pos] = False
            return True
        return False

    def refine(self, init: Optional[Sequence[bool]] = None) -> np.ndarray:
        sel = (np.ones(len(self.idx), dtype=bool) if init is None
               else np.array(init, dtype=bool))
        n = len(self.idx)
        seen = set()
        # Residues re-added by the step-7 gap fill are never removed again.
        # Without this, 1DGQ/1BQ0 oscillate between removing a residue (which
        # opens a small gap) and refilling it; the reference ends with the gap
        # filled.
        refilled = np.zeros(n, dtype=bool)
        while True:
            # Cycle guard (defensive): a repeated selection means no progress.
            key = sel.tobytes()
            if key in seen:
                break
            seen.add(key)
            r, _ = self._rmsd(sel)
            # -u RMSD threshold: once the current range is tight enough, stop
            # trimming. The reference implementation keeps removing mobile residues only while the
            # domain RMSD stays above `rmsd_threshold` (e.g. 1d3z collapses to
            # 0.096 A < 0.1 after its tail is removed, so trimming stops early).
            if r <= self.params.rmsd_threshold:
                break
            # step 2: remove isolated residues (two neighbours not selected)
            changed = False
            for pos in range(n):
                if not sel[pos]:
                    continue
                left = sel[pos - 1] if pos - 1 >= 0 else False
                right = sel[pos + 1] if pos + 1 < n else False
                if not left and not right:
                    sel[pos] = False
                    changed = True
            if changed:
                continue
            # step 3-5
            disp = self._displacements(sel)
            best_no_gap, best_gap = None, None
            for pos in range(n):
                if not sel[pos] or refilled[pos]:
                    continue
                if self._opens_gap(sel, pos):
                    if best_gap is None or disp[self.idx[pos]] > disp[self.idx[best_gap]]:
                        best_gap = pos
                else:
                    if best_no_gap is None or disp[self.idx[pos]] > disp[self.idx[best_no_gap]]:
                        best_no_gap = pos
            cands = []
            for pos in (best_no_gap, best_gap):
                if pos is None:
                    continue
                cands.append((self._decrease(sel, pos, r), pos))
            removed = False
            if cands:
                dr, pos = max(cands, key=lambda t: t[0])
                removed = self._condition_remove(sel, pos, r, dr)
            if not removed:
                # step 6: retry with the overall largest gap-weighted decrease
                best_pos, best_dr = None, -1.0
                for pos in range(n):
                    if not sel[pos] or refilled[pos]:
                        continue
                    dr = self._decrease(sel, pos, r)
                    if dr > best_dr:
                        best_dr, best_pos = dr, pos
                if best_pos is not None:
                    removed = self._condition_remove(sel, best_pos, r, best_dr)
            if removed:
                continue
            # step 7: fill small gaps (< g residues)
            before = sel.copy()
            if not self._fill_small_gaps(sel):
                break
            refilled |= sel & ~before
        return sel

    def _fill_small_gaps(self, sel: np.ndarray) -> bool:
        p = self.params
        changed = False
        i = 0
        n = len(self.idx)
        while i < n:
            if not sel[i]:
                j = i
                while j < n and not sel[j]:
                    j += 1
                width = j - i
                if 0 < i and j < n and width < p.gap_width:
                    for k in range(i, j):
                        sel[k] = True
                    changed = True
                i = j
            else:
                i += 1
        return changed


def _note_few_core(structure: Structure, n_core: int, params: Params) -> None:
    structure.few_core = True
    structure.warnings.append(
        f"*** WARNING: The number of core atoms ({n_core:3d}) was below the domain "
        f"threshold value ({params.min_cluster:3d}).")


def identify_domains(structure: Structure, params: Params) -> List[Domain]:
    """Stage A + B driver returning the domains with refined ranges.

    First-cut: a single domain spanning the full sequence (single-domain protein).
    Full Stage A (order-parameter core atoms + distance-variance Ward clustering)
    is layered in as `_stage_a_ranges`.
    """
    n = len(structure.residues)
    full = list(range(n))
    if params.method == "default":
        starts = _default_stage_a(structure, params)
    else:
        starts = [(rng, None) for rng in _stage_a_ranges(structure, params, full)]
    domains = []
    for rng, init in starts:
        d = evaluate(structure, params, rng, init)
        # Default method: with too few core atoms the whole-chain range is reported down to
        # 3 residues, as the reference does (a 3-residue range for a 6-residue chain).
        if d.n_residues >= params.min_cluster or (structure.few_core and d.n_residues >= 3):
            domains.append(d)
    if params.method == "default":
        if DEFAULT_DOMAIN_RATIO is not None and domains:
            # A domain far less ordered than the best-defined one is disorder, not a
            # domain (pieces of a disordered chain next to an ordered segment).
            best = min(d.rmsd for d in domains)
            domains = [d for d in domains if d.rmsd <= DEFAULT_DOMAIN_FLOOR
                       or d.rmsd <= DEFAULT_DOMAIN_RATIO * best]
        domains.sort(key=lambda d: (d.chains[0] if d.chains else "", d.ranges[0]))
    return domains


def identify_domains_per_state(structure: Structure, params: Params,
                               states) -> List[Tuple[object, List[Domain]]]:
    """The whole Stage A + B analysis run on each conformational state alone.

    Returns ``[(state, domains), ...]`` in the state order of `states`. A state
    is an ensemble in its own right, so its answer is a fresh run on just those
    conformers, not a filter of what the full bundle reported -- and the two
    genuinely differ. A full bundle reports an RMSD inflated by the motion
    *between* states, and the gap it leaves between two domains is often that
    inter-state hinge rather than disorder inside a state, so within a state the
    gap can refill and the domains merge.

    A state of a single conformer has no internal spread to measure, so it is
    left out.
    """
    out = []
    for state in states.states:
        if len(state.indices) < 2:
            continue
        out.append((state, identify_domains(
            subset_conformers(structure, state.indices), params)))
    return out


# Stage-A order-parameter cutoff. Rather than a fixed threshold, the paper selects
# the S_i (over all rotatable torsions in the structure) maximizing
#   Q_i = ((S_i - S_min)/(S_max - S_min)) * ((s - r_i + 1)/s),
# where r_i is the ascending rank of S_i (1 = smallest) and s the total number of
# torsions. This lets a per-structure cutoff drop genuinely mobile torsions (the
# noise clusters feeding Stage A) while retaining the rigid backbone.
STAGE_A_RIGID_FRAC = 0.07   # cluster is 'rigid' if intra-group variance <= this fraction
STAGE_A_SEP = 0.05          # minimal intra-variance increase to accept a split
STAGE_A_PAINT = False       # paint overlapping buffered ranges (reference rule; off: costs two_domain, 1aaf)
STAGE_A_RMSD_SPLIT = False  # domains from the RMSD split test (_rmsd_split_groups)
# Default method, Stage-B removal gate. The paper's text requires both decrease
# conditions (AND); the reference program accepts either (OR). With the default
# Stage A on the 134 evaluation cases, AND gives looser ranges (coverage 0.890,
# mean RMSD 0.708 A) while OR matches the published behaviour (coverage 0.822,
# RMSD 0.560 A vs the reference's 0.828 / 0.565), so OR is the default.
DEFAULT_AND_GATE = False
# Default method, torsions for the order parameters: "bb" (phi, psi), "bb_chi1",
# "chi" (phi, psi, all chi), "chi_h" (chi plus the _H_ROTORS hydrogen rotors). On the
# 134 evaluation cases all four agree (domain count 130-131, Jaccard 0.963-0.965), as the
# paper states; "chi" is the default because it does not depend on whether, or how, the
# input names hydrogens.
DEFAULT_TORSIONS = "chi"
# Split test (default method): split iff r(parent) - r(left) - r(right) >= SPLIT_TAU (A)
# and >= SPLIT_REL * the accepted parent split's value. SPLIT_TAU is the centre of the
# plateau of scripts/stagea/sens_split.py (reference domain count on 128-131/134 real
# cases over tau 0.06-0.20). SPLIT_REL was 0.6 until the improvement loop of
# scripts/heldout/PROTOCOL.md: it blocked clear sub-splits under a large parent split
# (cdk_c4d3: each chain r 1.7 A, halves 0.4 A, kept whole) and merged bodies on the
# held-out sets; 0 removes the rule.
SPLIT_TAU = 0.12
SPLIT_REL = 0.0
# Tree for the split test: "ward" = Ward linkage on V treated as distances;
# "average_sd" = average linkage on sqrt(V) (Angstrom). Both give the same domain counts
# except one case (2DHS: average_sd keeps one domain; 130 vs 131/134), so "ward".
SPLIT_LINKAGE = "ward"
# A node whose own split fails the test is still divided when an acceptable split
# exists below it: a child made of several bodies inflates r(child), so m is small
# although the node is not one rigid body (truth benchmark tail16_2dom). Splits found
# that way must leave parts of at least ceil(C/8) core atoms (the paper's granularity:
# average domain size above one eighth of the core atoms), so short rigid elements
# inside a larger domain do not become domains (1a03 helices). Stage-A domain counts:
# real 130/134, synthetic 68/72, ground truth 13/13 (without: 131, 68, 12).
SPLIT_DEEP = True
SPLIT_DEEP_MINFRAC = 0.125
# A node too small to split into two domains (both halves under mu atoms) is discarded
# when it is not one rigid body. "sum": r(node) - r(left) - r(right) >= SPLIT_TAU, i.e. it
# is two bodies; that misses a uniformly flexible segment, whose halves are flexible too
# (disordered tails kept as domains on the held-out benchmark, set 1). "max":
# r(node) - max(r(left), r(right)) >= SPLIT_TAU: even the more disordered half is more
# ordered than the whole by the RMSD that defines a domain split. "max" discards whatever
# "sum" discards. Disordered tail pieces lose 0.5-0.8 A this way, small rigid domains of
# real ensembles at most 0.08 A. (A scale-free ratio, r(left) + r(right) < 1.5 r(node),
# fixed the benchmark but discarded small rigid domains of real ensembles, whose halves
# are also more ordered than the whole, at a level of a few hundredths of an Angstrom.)
FRAG_TEST = "max"
# Default start ranges: "buffer" = each domain's core span extended by m residues;
# "partition" = every residue of the chain goes to the domain it moves with (its core
# domain, else the one in whose frame it moves least), and Stage B trims.
DEFAULT_SPAN = "buffer"
# Largest RMSD (A) of a group accepted as one rigid domain; a larger unsplit group is
# searched for rigid parts (None: no limit).
SPLIT_RMAX: Optional[float] = None
# Reported domains with an RMSD above DEFAULT_DOMAIN_FLOOR (A) and above DEFAULT_DOMAIN_RATIO
# times the lowest domain RMSD of the ensemble are dropped (None: keep all).
DEFAULT_DOMAIN_RATIO: Optional[float] = None
DEFAULT_DOMAIN_FLOOR = 1.0
# After the refinement, residues of the start range that move with the domain (mean
# displacement in its frame at most DEFAULT_READMIT x that of its residues) are selected
# again (None: off). Counters over-trimming of noisy but rigid residues.
DEFAULT_READMIT: Optional[float] = None
# Second split criterion (None: off): the node's Ward merge height is SPLIT_GAP times
# those of its children (two blocks in V) and the two-body fit explains an RMS motion of
# at least SPLIT_GAP_D A. Detects small hinge motions whose RMSD gain m is hidden by
# the noise of the two parts.
# Chosen on the development data (held-out sets 1-2, real ensembles): hinge-sweep
# detections 29/31 -> 39/41 of 56, real-ensemble results unchanged.
SPLIT_GAP: Optional[float] = 7.0
SPLIT_GAP_D = 0.5
# When every group is discarded as flexible, keep the discarded fragment with the lowest
# RMSD instead of falling back to the whole chain (1Y4E: 0.33 A instead of 2.19 A).
FRAG_FALLBACK = True


def _s_cut_q(structure: Structure) -> float:
    """Q_i-optimal order-parameter cutoff for the whole structure."""
    o = order_parameters(structure)
    vals = np.array([v for d in o for v in d.values()])
    s = vals.size
    if s == 0:
        return 0.0
    smin, smax = vals.min(), vals.max()
    if smax - smin < 1e-12:
        return smin
    order = np.argsort(vals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(s) + 1
    q = ((vals - smin) / (smax - smin)) * ((s - ranks + 1) / s)
    return float(vals[int(np.argmax(q))])


def _core_residue_indices(structure: Structure, s_cut: float) -> List[int]:
    """0-based residue indices whose CA is a core atom (some torsion S > s_cut)."""
    o = order_parameters(structure)
    return [i for i, d in enumerate(o) if any(v > s_cut for v in d.values())]


def _ca_variance_matrix(structure: Structure, core_idx: List[int]) -> np.ndarray:
    """V_ij = variance over conformers of the CA-CA distance between core atoms."""
    conformers = structure.conformers
    residues = structure.residues
    k = len(core_idx)
    ca = []
    for i in core_idx:
        r = residues[i]
        ca.append(np.array([c[(r.chain, r.num)]["CA"] for c in conformers]))
    V = np.zeros((k, k))
    for a in range(k):
        for b in range(k):
            V[a, b] = float(np.linalg.norm(ca[a] - ca[b], axis=1).var())
    return V


def _split_rigid_domains(structure: Structure, core_idx: List[int], mu: int) \
        -> List[Tuple[int, int]]:
    """(lo, hi) residue-index span of each group from `_rigid_domain_groups`."""
    ranges = [(min(g), max(g)) for g in _rigid_domain_groups(structure, core_idx, mu)]
    ranges.sort()
    return ranges


def _rigid_domain_groups(structure: Structure, core_idx: List[int], mu: int) \
        -> List[List[int]]:
    """Recursively split core atoms into well-separated rigid domains.

    A Ward-hierarchical cut into two groups is accepted as a domain split only if
    each resulting group is itself a rigid body (intra-group variance below a
    fraction of the total) and at least `mu` elements. A top-level cut that leaves
    exactly one large non-rigid group alongside sub-mu noise (e.g. a mobile
    terminus cluster) is still pursued: the noise is dropped and the remainder is
    re-split, so a single protein whose rigid bodies sit inside one flexible
    cluster can be partitioned (1aaf). Returns the domains as lists of core
    residue indices.
    """
    def recurse(idx: np.ndarray) -> List[List[int]]:
        n = len(idx)
        if n < mu:
            return []                       # sub-mu cluster -> noise, drop
        if n < 2 * mu:
            return [list(idx)]              # too small to split further
        sub = V[np.ix_(idx, idx)]
        if sub.sum() / 2.0 <= 1e-12:
            return [list(idx)]
        Z = linkage(squareform(sub, checks=False), method="ward")
        labels = fcluster(Z, 2, criterion="maxclust")
        sep = float(Z[:, 2][-1])
        total_sub = sub.sum() / 2.0
        children = []                       # (local_indices, is_rigid)
        for c in np.unique(labels):
            g = np.where(labels == c)[0]
            intra = sub[np.ix_(g, g)].sum() / 2.0 / total_sub
            children.append((g, intra <= STAGE_A_RIGID_FRAC))
        clean = sep > STAGE_A_SEP and all(
            rigid or len(g) < mu for (g, rigid) in children)
        n_rigid_big = sum(1 for (g, rigid) in children if rigid and len(g) >= mu)
        if clean and n_rigid_big >= 2:
            out = []
            for (g, rigid) in children:
                out.extend(recurse(idx[g]))
            return out
        # Not clean: pursue a split that dropped noise. If exactly one large
        # non-rigid child remains and every other child is sub-mu noise, drop the
        # noise and recurse into the large remainder (the 1aaf pattern).
        big_not_rigid = [g for (g, rigid) in children if len(g) >= mu and not rigid]
        rest_noise = all(len(g) < mu for (g, rigid) in children
                         if not (len(g) >= mu and not rigid))
        if len(big_not_rigid) == 1 and rest_noise:
            g = big_not_rigid[0]
            out = recurse(idx[g])
            if out:
                return out
        return [list(idx)]

    k = len(core_idx)
    if k < 2 * mu:
        return [list(core_idx)]
    V = _ca_variance_matrix(structure, core_idx)
    if V.sum() / 2.0 <= 1e-12:
        return [list(core_idx)]
    return [[core_idx[i] for i in g] for g in recurse(np.arange(k))]


def _rmsd_split_groups(structure: Structure, core_idx: List[int], mu: int,
                       info: Optional[Dict[str, bool]] = None) -> List[List[int]]:
    """Domains by a top-down RMSD split test on the Ward tree of the core atoms.

    m = r(parent) - r(left) - r(right), r = RMSD to the mean of the cluster's N/CA/C
    atoms (single fit to conformer 1). Derived on synthetic hinge ensembles, AISAR
    bundles and the real readouts (scripts/stagea/probe_rootsplit.py, probe_absdelta.py;
    threshold sensitivity in scripts/stagea/sens_split.py):
      * both children >= mu: split iff m >= SPLIT_TAU and m >= SPLIT_REL * m of the
        accepted parent split; otherwise the node is one domain;
      * one child < mu: the small child is discarded as noise and the big one examined
        (a sub-mu fragment never decides the split of the rest);
      * both children < mu: the node is discarded if halving it lowers the RMSD even of
        its more disordered half by SPLIT_TAU, r(node) - max(r(left), r(right)) >=
        SPLIT_TAU (internally flexible segments; FRAG_TEST "sum": m >= SPLIT_TAU),
        otherwise it is one domain.
    """
    k = len(core_idx)
    if k < 2:
        return [list(core_idx)]
    V = _ca_variance_matrix(structure, core_idx)
    if V.sum() / 2.0 <= 1e-12:
        return [list(core_idx)]
    if SPLIT_LINKAGE == "average_sd":
        D = np.sqrt(np.maximum(V, 0.0))      # distance standard deviation (A)
        np.fill_diagonal(D, 0.0)
        Z = linkage(squareform((D + D.T) / 2, checks=False), method="average")
    else:
        Z = linkage(squareform(V, checks=False), method="ward")
    root = to_tree(Z)
    res = structure.residues
    X = np.array([[[c[(res[i].chain, res[i].num)][a] for a in ("N", "CA", "C")]
                   for i in core_idx] for c in structure.conformers])
    cache: Dict[Tuple[int, ...], float] = {}

    def r(g: Tuple[int, ...]) -> float:
        if g not in cache:
            ens = X[:, list(g)].reshape(X.shape[0], -1, 3)
            cache[g] = rmsd_to_mean(ens, ref_seq=ens[0], iters=1)
        return cache[g]

    out: List[Tuple[int, ...]] = []
    dropped: List[Tuple[int, ...]] = []     # fragments discarded as flexible

    def accepted(m: float, parent_m: Optional[float]) -> bool:
        return m >= SPLIT_TAU and (parent_m is None or m >= SPLIT_REL * parent_m)

    def blocky(node, ids, L, R) -> bool:
        """Two-block evidence: the Ward merge height of the node is SPLIT_GAP times
        those of its children and the two-body model explains an RMS motion
        sqrt(r^2 - r2^2) >= SPLIT_GAP_D (r2 the pooled residual of the two parts)."""
        if SPLIT_GAP is None:
            return False
        gap = node.dist / max(node.left.dist, node.right.dist, 1e-12)
        if gap < SPLIT_GAP:
            return False
        r2 = (len(L) * r(L) ** 2 + len(R) * r(R) ** 2) / len(ids)
        return r(ids) ** 2 - r2 >= SPLIT_GAP_D ** 2

    # Splits found below a failed node must leave parts of at least deep_min atoms
    # (SPLIT_DEEP_MINFRAC of the core, at least mu): small rigid elements inside a
    # larger domain are not domains of their own.
    deep_min = max(mu, int(np.ceil(SPLIT_DEEP_MINFRAC * k)))

    def has_split(node, parent_m: Optional[float]) -> bool:
        """Whether an acceptable domain split exists at `node` or below it."""
        ids = tuple(sorted(node.pre_order()))
        if node.is_leaf() or len(ids) < mu:
            return False
        L = tuple(sorted(node.left.pre_order()))
        R = tuple(sorted(node.right.pre_order()))
        if len(L) >= deep_min and len(R) >= deep_min and \
                (accepted(r(ids) - r(L) - r(R), parent_m) or blocky(node, ids, L, R)):
            return True
        big = [c for c, g in ((node.left, L), (node.right, R)) if len(g) >= mu]
        if len(big) < 2:
            return bool(big) and has_split(big[0], parent_m)
        return has_split(node.left, parent_m) or has_split(node.right, parent_m)

    def walk(node, parent_m: Optional[float]) -> None:
        ids = tuple(sorted(node.pre_order()))
        if node.is_leaf() or len(ids) < mu:
            if len(ids) >= mu:
                out.append(ids)
            return
        L = tuple(sorted(node.left.pre_order()))
        R = tuple(sorted(node.right.pre_order()))
        big = [c for c, g in ((node.left, L), (node.right, R)) if len(g) >= mu]
        if len(big) == 1:
            walk(big[0], parent_m)
            return
        m = r(ids) - r(L) - r(R)
        if len(big) == 2:
            if accepted(m, parent_m) or blocky(node, ids, L, R):
                walk(node.left, m)
                walk(node.right, m)
            elif SPLIT_DEEP and (has_split(node.left, parent_m) or
                                 has_split(node.right, parent_m)):
                # A child that is itself several bodies inflates r(child), so m can
                # be small although the node is not rigid: divide it anyway.
                walk(node.left, parent_m)
                walk(node.right, parent_m)
            elif SPLIT_RMAX is not None and r(ids) > SPLIT_RMAX:
                # Too disordered to be one rigid body (a disordered chain whose torsions
                # all pass the relative order cut): look for rigid parts inside.
                walk(node.left, parent_m)
                walk(node.right, parent_m)
            else:
                out.append(ids)
        elif FRAG_TEST == "sum":
            if m < SPLIT_TAU:
                out.append(ids)
        elif r(ids) - max(r(L), r(R)) < SPLIT_TAU and \
                (SPLIT_RMAX is None or r(ids) <= SPLIT_RMAX):
            out.append(ids)
        else:
            dropped.append(ids)

    walk(root, None)
    if not out and dropped and FRAG_FALLBACK:
        # Nothing passed as rigid: the most ordered discarded fragment is the best
        # superposition core there is (1Y4E), better than the whole flexible chain.
        best = min(dropped, key=r)
        if SPLIT_RMAX is None or r(best) <= SPLIT_RMAX:
            out.append(best)
            if info is not None:
                info["fallback"] = True
    return [[core_idx[i] for i in g] for g in out]


def _paint_domains(structure: Structure, groups: List[List[int]], params: Params) \
        -> List[List[int]]:
    """Resolve overlapping buffered domain ranges (synthetic-ensemble probes).

    Each domain's buffered range [core_min - b, core_max + b] is painted onto the
    chain in ascending order of the group's mean intra-cluster CA distance
    variance, later paint overwriting earlier, so the looser body wins every
    overlap. A domain left with fewer than `min_cluster` residues is dropped.
    """
    n = len(structure.residues)
    keyed = []
    for g in groups:
        m = len(g)
        V = _ca_variance_matrix(structure, g)
        keyed.append((V.sum() / (m * (m - 1)) if m > 1 else 0.0,
                      max(0, min(g) - params.buffer), min(n - 1, max(g) + params.buffer)))
    owner = [-1] * n
    for d, (_, lo, hi) in enumerate(sorted(keyed)):
        for i in range(lo, hi + 1):
            owner[i] = d
    out = []
    for d in range(len(keyed)):
        res = [i for i in range(n) if owner[i] == d]
        if len(res) >= params.min_cluster:
            out.append(res)
    out.sort()
    return out


def _stage_a_ranges(structure: Structure, params: Params, full: Sequence[int]):
    """Stage A domain-boundary ranges (list of residue-index lists).

    Computes order parameters -> core CA atoms -> distance-variance matrix ->
    Ward clustering -> rigid-domain split (see `_split_rigid_domains`), extending
    each domain by the boundary buffer. Single-domain proteins return the whole
    sequence as one range.
    """
    n = len(structure.residues)
    core = _core_residue_indices(structure, _s_cut_q(structure))
    if len(core) < params.min_cluster:
        # No few-core warning here: this core count differs from the reference's on
        # some ensembles (1DLZ: 7 here, at least 8 there), and reference mode comes first.
        return [list(full)]
    if STAGE_A_RMSD_SPLIT:
        groups = _rmsd_split_groups(structure, core, params.min_cluster)
    else:
        groups = _rigid_domain_groups(structure, core, params.min_cluster)
    if len(groups) <= 1:
        return [list(full)]
    if STAGE_A_PAINT:
        return _paint_domains(structure, groups, params) or [list(full)]
    out = []
    for lo, hi in sorted((min(g), max(g)) for g in groups):
        out.append(list(range(max(0, lo - params.buffer),
                              min(n, hi + params.buffer + 1))))
    return out


def _published_s_cut(order: List[Dict[str, float]]) -> float:
    """S_cut maximising the published Q_i = (s-1)(S_i-S_min)/(S_max-S_min) - r_i.

    r_i is the ascending rank of S_i (1 = smallest) and s the number of torsions;
    ties of the maximum go to the highest S. On three-level synthetic populations
    this reproduces the reference cut exactly (scripts/synthetic/levelprobe.py).
    """
    vals = np.sort(np.array([v for d in order for v in d.values()]))
    s = vals.size
    if s == 0:
        return 0.0
    lo, hi = vals[0], vals[-1]
    if hi - lo < 1e-12:
        return float(lo)
    q = (s - 1) * (vals - lo) / (hi - lo) - np.arange(1, s + 1)
    return float(vals[np.where(q >= q.max() - 1e-12)[0].max()])


def _frame_displacements(structure: Structure, core_idx: List[int]) -> np.ndarray:
    """Per-residue backbone displacement after superposing on `core_idx`.

    Every conformer is fitted (N, CA, C of the core residues) onto conformer 1 and
    then onto the mean of that fit; the result for a residue is the average over
    conformers and its N, CA, C atoms of the distance to the mean position.
    Residues lacking a backbone atom get +inf.
    """
    res = structure.residues
    confs = structure.conformers
    X = np.full((len(confs), len(res), 3, 3), np.nan)
    for k, c in enumerate(confs):
        for i, r in enumerate(res):
            atoms = c.get((r.chain, r.num), {})
            for j, a in enumerate(BACKBONE):
                if a in atoms:
                    X[k, i, j] = atoms[a]
    fit = X[:, core_idx].reshape(len(confs), -1, 3)
    ref = fit[0]
    for _ in range(2):
        moved = []
        for k in range(len(confs)):
            R = _kabsch_matrix(fit[k], ref)
            moved.append((X[k] - fit[k].mean(axis=0)) @ R.T + ref.mean(axis=0))
        moved = np.array(moved)
        ref = moved[:, core_idx].reshape(len(confs), -1, 3).mean(axis=0)
    mean = moved.mean(axis=0)
    d = np.linalg.norm(moved - mean, axis=3).mean(axis=(0, 2))
    return np.where(np.isnan(d), np.inf, d)


def _kabsch_matrix(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rotation matrix superposing point set P onto Q (same order)."""
    H = (P - P.mean(axis=0)).T @ (Q - Q.mean(axis=0))
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _default_domain_ranges(structure: Structure, groups: List[List[int]], params: Params,
                         span: str = "buffer") -> List[Tuple[List[int], List[bool]]]:
    """Start ranges for Stage B, built the same way for one or several domains.

    Each domain claims its core span extended by m residues on both sides (span
    "partition": every residue of the chain is claimed by every domain). A residue
    claimed by several domains goes to the domain whose core atom it is, otherwise
    to the domain in whose frame it moves least (`_frame_displacements`). Returns
    (contiguous span, initial selection) per domain.
    """
    n = len(structure.residues)
    m = params.buffer
    spans = [(max(0, min(g) - m), min(n - 1, max(g) + m)) for g in groups]
    claims: Dict[int, List[int]] = {}
    if span == "partition":
        # Every residue is a candidate for the domain it moves with; Stage B decides.
        claims = {i: list(range(len(groups))) for i in range(n)}
    for d, (lo, hi) in enumerate(spans):
        for i in range(lo, hi + 1):
            if i not in claims or d not in claims[i]:
                claims.setdefault(i, []).append(d)
    core_owner = {i: d for d, g in enumerate(groups) for i in g}
    disp: Dict[int, np.ndarray] = {}
    owner: Dict[int, int] = {}
    for i, ds in claims.items():
        if len(ds) == 1:
            owner[i] = ds[0]
        elif core_owner.get(i) in ds:
            owner[i] = core_owner[i]
        else:
            for d in ds:
                if d not in disp:
                    disp[d] = _frame_displacements(structure, groups[d])
            owner[i] = min(ds, key=lambda d: (disp[d][i], d))
    out = []
    for d, (lo, hi) in enumerate(spans):
        if span == "partition":
            mine = [i for i in range(n) if owner[i] == d]
            lo, hi = min(mine + [lo]), max(mine + [hi])
        rng = list(range(lo, hi + 1))
        out.append((rng, [owner[i] == d for i in rng]))
    return out


def _default_stage_a(structure: Structure, params: Params)\
        -> List[Tuple[List[int], Optional[List[bool]]]]:
    """RESUPER Stage A (params.method == "default").

    Core residues: any own torsion with S >= the published Q_i cutoff. Domains: the
    RMSD-decrease split test on the Ward tree of the core distance variances
    (`_rmsd_split_groups`). Ranges: `_default_domain_ranges`, also for a single
    domain. Too few core atoms -> one domain over the whole chain.
    """
    n = len(structure.residues)
    order = order_parameters(structure, h_rotors=DEFAULT_TORSIONS == "chi_h")
    keep = {"bb": ("phi", "psi"), "bb_chi1": ("phi", "psi", "chi1")}.get(DEFAULT_TORSIONS)
    if keep:
        order = [{k: v for k, v in d.items() if k in keep} for d in order]
    cut = _published_s_cut(order)
    core = [i for i, d in enumerate(order) if any(v >= cut for v in d.values())]
    if len(core) < params.min_cluster:
        _note_few_core(structure, len(core), params)
        return [(list(range(n)), None)]
    info: Dict[str, bool] = {}
    groups = _rmsd_split_groups(structure, core, params.min_cluster, info)
    if not groups:
        return [(list(range(n)), None)]
    return _default_domain_ranges(structure, groups, params,
                                span="buffer" if info.get("fallback") else DEFAULT_SPAN)


def _readmit(structure: Structure, domain_range: Sequence[int], sel: Sequence[bool],
             init: Optional[Sequence[bool]]) -> List[bool]:
    """Give back residues the refinement removed but that move with the domain: in the
    frame of the selected residues, a residue of the start range (owned by the domain)
    whose mean backbone displacement is at most DEFAULT_READMIT times the mean
    displacement of the selected residues is selected again (repeated until stable)."""
    sel = list(sel)
    for _ in range(5):
        chosen = [domain_range[p] for p in range(len(sel)) if sel[p]]
        if len(chosen) < 3:
            return sel
        disp = _frame_displacements(structure, chosen)
        limit = DEFAULT_READMIT * float(np.mean(disp[chosen]))
        changed = False
        for p, i in enumerate(domain_range):
            if not sel[p] and (init is None or init[p]) and disp[i] <= limit:
                sel[p] = True
                changed = True
        if not changed:
            break
    return sel


def evaluate(structure: Structure, params: Params, domain_range: Sequence[int],
             init: Optional[Sequence[bool]] = None) -> Domain:
    refiner = RangeRefiner(structure, params, domain_range)
    sel = refiner.refine(init)
    if params.method == "default" and DEFAULT_READMIT is not None:
        sel = _readmit(structure, domain_range, sel, init)
    n = len(sel)
    sel_abs = np.zeros(len(structure.residues), dtype=bool)
    for pos in range(n):
        if sel[pos]:
            sel_abs[domain_range[pos]] = True
    segs = _chain_segments(sel_abs, structure.residues)
    ranges = [(lo, hi) for _, lo, hi in segs]
    gap_count = len(segs) - 1 if segs else 0
    # Multi-chain structures label every segment with its chain ID (e.g. A4..A92).
    multichain = len({r.chain for r in structure.residues}) > 1
    chains = [c if multichain else "" for c, _, _ in segs]
    n_res = int(sel_abs.sum())
    r = refiner._reported_rmsd(sel)
    return Domain(ranges=ranges, rmsd=r, n_residues=n_res, n_gaps=gap_count,
                  chains=chains)
