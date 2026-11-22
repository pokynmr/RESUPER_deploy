"""Command-line interface for RESUPER.

    resuper [PARAMETERS] <input file>

RESUPER identifies ordered regions (domains) and the residues to superimpose in
a multi-conformer protein ensemble. The default method is RESUPER's own
interpretable domain/range procedure. `--method reference` reproduces the
report of the CYRANGE 2.0 reference implementation so that the two methods can
be compared side by side on the same input.

Parameters (defaults): -a absdecrease 1.6, -b buffer 3, -d decrease 1.2,
-g gap 0.4, -u RMSD threshold 0.1, -gw gap width 3, -l library, -m min cluster 8,
-r residue range.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional, Sequence

from .core import Params, Domain, identify_domains, identify_domains_per_state
from .pdb import parse_pdb
from .states import detect_states
from .writer import superpose_ensemble, superpose_to_mean, write_cif, write_pdb

STAR = "*" * 44
FILLER = " "


def _range_lines(domains: List[Domain], min_cluster: int = 8) -> List[str]:
    """The `Optimal range ...` lines for one set of domains.

    Shared by the ensemble report and the per-state report so both print a
    range the same way."""
    lines = []
    if domains:
        for d in domains:
            chains = d.chains or [""] * len(d.ranges)
            ranges = ", ".join(f"{c}{a}" if a == b else f"{c}{a}..{c}{b}"
                               for c, (a, b) in zip(chains, d.ranges))
            gword = "gap" if d.n_gaps == 1 else "gaps"
            rword = "residue" if d.n_residues == 1 else "residues"
            lines.append(
                f"    Optimal range {ranges}: RMSD {d.rmsd:.2f} A, "
                f"{d.n_gaps} {gword}, {d.n_residues} {rword}.")
    else:
        lines.append(
            f"*** WARNING: Only domains with <{min_cluster:3d} residues could be identified.")
    return lines


def format_report(pdb_name: str, n_conformers: int, domains: List[Domain],
                  min_cluster: int = 8, warnings: Sequence[str] = (),
                  reference: bool = False) -> str:
    if reference:
        header = [STAR, " " * 14 + "CYRANGE (v. 2.0)", FILLER,
                  " " * 4 + "Copyright (c) 2002-13 Peter Guntert.", STAR]
    else:
        header = [STAR, " " * 14 + "RESUPER (v. 1.0)", FILLER,
                  " " * 4 + "Ordered regions and superposition ranges.", STAR]
    lines = header + [
        f"    PDB coordinate file \"{pdb_name}\" read, {n_conformers} rigid conformers.",
        FILLER] + list(warnings)
    lines += _range_lines(domains, min_cluster)
    lines.append(FILLER)
    return "\n".join(lines) + "\n"


def _conformer_runs(members) -> str:
    """Compact conformer numbers into runs: 1, 3, 6, 7, 10 -> "1, 3, 6-7, 10".

    Only genuinely consecutive numbers are collapsed.  A two-state ensemble
    interleaves its conformers (7SA5 state 1 is 1, 3, 6, 7, ...), so showing a
    single first-last span would print overlapping ranges and hide the split."""
    nums = sorted(int(m) for m in members)
    parts = []
    start = prev = nums[0]
    for m in nums[1:] + [None]:
        if m != prev + 1:
            parts.append(f"{start}-{prev}" if prev > start else str(start))
            start = m
        prev = m
    return ", ".join(parts)


def format_states_report(n_conformers: int, report) -> str:
    """Human-readable multiconformer detection line appended to the report."""
    populations = ", ".join(f"{s.population:.2f}" for s in report.states)
    parts = [
        f"conformers {_conformer_runs(s.members)} -> state {i + 1}"
        for i, s in enumerate(report.states)
    ]
    groups = ", ".join(parts)
    return (
        f"    Conformational states: {report.n_states} "
        f"(populations {populations}); {groups}."
    )


def format_per_state_report(per_state, min_cluster: int = 8) -> str:
    """Per-state ordered ranges, printed under the whole-ensemble report.

    `per_state` is `core.identify_domains_per_state`'s ``[(state, domains)]``.
    The ensemble-wide lines above stay in the report: a state's ranges answer a
    different question (what is ordered *within* that state), and only the pair
    shows whether an ensemble-wide gap is inter-state motion or disorder.
    """
    lines = []
    for i, (state, domains) in enumerate(per_state, start=1):
        lines.append(f"    State {i} (conformers {_conformer_runs(state.members)}):")
        lines += _range_lines(domains, min_cluster)
    return "\n".join(lines) + "\n"


def _coord_format(path: Optional[str]) -> str:
    if path:
        ext = os.path.splitext(path)[1].lower()
        if ext in (".cif", ".mmcif"):
            return "cif"
    return "pdb"


def _write_outputs(structure, ns, states) -> None:
    """Write a superposed ensemble (and per-state splits) onto the chosen
    reference (first model, or the iterative ensemble mean)."""
    if ns.superpose_ref == "mean":
        superposed = superpose_to_mean(structure)
    else:
        superposed = superpose_ensemble(structure, 0)
    stem = os.path.splitext(os.path.basename(ns.pdb))[0]
    ens_paths = ns.write_ensemble or []
    # Split format follows a lone .cif ensemble file, else PDB.
    split_fmt = _coord_format(ens_paths[0]) if len(ens_paths) == 1 else "pdb"
    ext = ".cif" if split_fmt == "cif" else ".pdb"
    write = write_cif if split_fmt == "cif" else write_pdb
    for path in ens_paths:
        (write_cif if _coord_format(path) == "cif" else write_pdb)(
            path, structure, superposed)
        print(f"*** wrote superposed ensemble to {path}", file=sys.stderr)
    if ns.split_states:
        if states is None:
            states = detect_states(structure)
        os.makedirs(ns.split_states, exist_ok=True)
        for i, st in enumerate(states.states, start=1):
            path = os.path.join(ns.split_states, f"{stem}_state{i}{ext}")
            write(path, structure, superposed, st.indices)
            print(f"*** wrote state {i} ({len(st.members)} conformers) to {path}",
                  file=sys.stderr)


def parse_args(argv: List[str]) -> tuple:
    p = argparse.ArgumentParser(add_help=False, prog="resuper")
    p.add_argument("-a", dest="absdecrease", type=float, default=1.6)
    p.add_argument("-b", dest="buffer", type=int, default=3)
    p.add_argument("-d", dest="decrease", type=float, default=1.2)
    p.add_argument("-g", dest="gap", type=float, default=0.4)
    p.add_argument("-u", dest="rmsd_threshold", type=float, default=0.1)
    p.add_argument("-gw", dest="gap_width", type=int, default=3)
    p.add_argument("-l", dest="lib")
    p.add_argument("-m", dest="min_cluster", type=int, default=8)
    p.add_argument("-r", dest="residue_range")
    p.add_argument("-h", dest="help", action="store_true")
    p.add_argument("--states", action="store_true",
                   help="detect and report the number of conformational "
                        "states (multiconformer analysis)")
    p.add_argument("--per-state", action="store_true",
                   help="analyse each conformational state as an ensemble of "
                        "its own and report its ranges (implies --states)")
    p.add_argument("--write-ensemble", metavar="FILE", action="append",
                   help="write the ensemble superposed onto the first model as "
                        "PDB (.pdb) or mmCIF (.cif) — all parsed atoms; "
                        "repeat to emit several files")
    p.add_argument("--split-states", metavar="DIR",
                   help="write one superposed PDB/CIF per detected "
                        "conformational state into DIR (requires --states)")
    p.add_argument("--superpose-ref", choices=("first", "mean"), default="first",
                   help="superpose the ensemble onto the first model (default) "
                        "or onto the iterative ensemble mean")
    p.add_argument("--method", choices=("default", "reference"), default="default",
                   help="default: RESUPER's own method (recommended); "
                        "reference: reproduce the CYRANGE 2.0 report for comparison")
    p.add_argument("pdb")
    ns = p.parse_args(argv)
    if ns.help:
        return None, ns
    # Per-state ranges only mean something once the states exist, and both
    # lines come from the same detection, so one flag implies the other rather
    # than leaving a "--per-state without --states" combination to report
    # nothing.
    if ns.per_state:
        ns.states = True
    params = Params(absdecrease=ns.absdecrease, buffer=ns.buffer, decrease=ns.decrease,
                    gap=ns.gap, rmsd_threshold=ns.rmsd_threshold,
                    gap_width=ns.gap_width, lib=ns.lib, min_cluster=ns.min_cluster,
                    residue_range=ns.residue_range, method=ns.method)
    return params, ns


HELP = r"""
  ****************** RESUPER HELP ******************

      Usage: resuper [PARAMETERS] <input file>

        The input file must be in PDB format.

   Optional Parameters (default values are given):
    -a: absdecrease (real; 1.6)
    -b: buffer (integer; 3)
    -d: decrease (real; 1.2)
    -g: gap (real; 0.4)
    -u: RMSD threshold (real; 0.1)
   -gw: gap width (integer; 3)
    -l: path to custom library file (internal lib)
    -m: minimum cluster size (integer; 8)
    -r: residue range (all residues)
    --method: default | reference
        default   RESUPER's own ordered-region and range method (recommended)
        reference reproduce the CYRANGE 2.0 report for side-by-side comparison
    --states  report the number of conformational states in the ensemble
    --per-state  report the ordered ranges of each conformational state
                 separately, as an ensemble of its own (implies --states)
    --write-ensemble FILE   write the superposed ensemble (PDB or mmCIF)

   EXAMPLE: resuper -g 0.3 -r 1-40,60-90 prot.pdb

  **************************************************
"""


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print("\n*** No input file was provided! ***\n    Usage: resuper [PARAMETERS] <input file>")
        print('    Use "resuper -h" for a help overview.')
        return 0
    if argv == ["-h"] or (argv and argv[0] == "-h"):
        print(HELP)
        return 0
    params, ns = parse_args(argv)
    try:
        structure = parse_pdb(ns.pdb)
        if not structure.residues:
            print("*** FATAL ERROR: No amino acid residues in sequence.")
            return 1
        if structure.n_models < 2:
            print("*** FATAL ERROR: Less than two structures selected.")
            return 1
        domains = identify_domains(structure, params)
    except Exception as e:      # never a traceback: report the error and fail
        print(f"*** FATAL ERROR: {type(e).__name__}: {e}")
        return 1
    if structure.dropped:
        print("*** WARNING: residues without N, CA and C in every conformer ignored: "
              + ", ".join(f"{c}{n}" for c, n in structure.dropped), file=sys.stderr)
    report = format_report(ns.pdb, structure.n_models, domains, params.min_cluster,
                           structure.warnings,
                           reference=params.method == "reference")
    states = None
    if ns.states:
        states = detect_states(structure)
        report += format_states_report(structure.n_models, states) + "\n"
    if ns.per_state:
        try:
            per_state = identify_domains_per_state(structure, params, states)
        except Exception as e:      # never a traceback: report the error and fail
            print(f"*** FATAL ERROR: {type(e).__name__}: {e}")
            return 1
        report += format_per_state_report(per_state, params.min_cluster)
    if ns.write_ensemble or ns.split_states:
        try:
            _write_outputs(structure, ns, states)
        except Exception as e:      # never a traceback: report and fail
            print(f"*** FATAL ERROR: {type(e).__name__}: {e}")
            return 1
    sys.stdout.write(report)
    return 0
