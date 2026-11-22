"""Confirmation sample: draw and freeze (decision S11).

Draws from the Phase 3 candidate pool (groups by rule S6), excluding the 60
Phase 3 candidates: A 45 at random, B all remaining, D all remaining, C 45 and
E 45 at random (numpy seed 20261113). Writes the list and its SHA-256 before any
structure is downloaded. Offline: uses only results/phase3/candidate_pool.csv.

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/confirm_select.py
Output: selection_study/results/confirm/candidates_confirm.csv, FROZEN_confirm.txt
"""
import csv
import hashlib
import os
import sys
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
OUT = os.path.join(S, "results", "confirm")
SEED = 20261113
SIZES = {"A": 45, "B": None, "D": None, "C": 45, "E": 45}   # None = all remaining


def main():
    if os.path.exists(os.path.join(OUT, "FROZEN_confirm.txt")) and "--redraw" not in sys.argv:
        print("Frozen selection exists (FROZEN_confirm.txt); not redrawing. The PDB changes over time, so a redraw "
              "would select different entries. Pass --redraw to draw a new selection.")
        return 0
    os.makedirs(OUT, exist_ok=True)
    pool = list(csv.DictReader(open(os.path.join(S, "results", "phase3", "candidate_pool.csv"))))
    used = {r["pdb_id"] for r in csv.DictReader(open(os.path.join(S, "results", "phase3",
                                                                  "candidates_selected.csv")))}
    rng = np.random.default_rng(SEED)
    chosen, report = [], []
    for g in "ABDCE":                                  # fixed order for the seeded draws
        members = sorted((r for r in pool if r["group"] == g and r["pdb_id"] not in used),
                         key=lambda r: r["pdb_id"])
        k = len(members) if SIZES[g] is None else min(SIZES[g], len(members))
        idx = sorted(rng.choice(len(members), size=k, replace=False)) if k else []
        chosen += [members[i] for i in idx]
        report.append(f"{g}: {k} of {len(members)} remaining")
    fields = ["group", "pdb_id", "protein", "organism", "length", "n_models", "chain", "group_reason",
              "n_dom_CATH", "n_dom_SCOP", "n_dom_ECOD", "residues_outside_domains", "helices",
              "n_strands", "title", "released", "citation", "pmid", "doi"]
    path = os.path.join(OUT, "candidates_confirm.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(chosen, key=lambda r: ("BADCE".index(r["group"]), r["pdb_id"])))
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(os.path.join(OUT, "FROZEN_confirm.txt"), "w") as fh:
        fh.write(f"candidates_confirm.csv sha256 {digest}\nfrozen {datetime.now(timezone.utc).isoformat()}\n"
                 f"seed {SEED}; {'; '.join(report)}; total {len(chosen)}\n")
    print("; ".join(report), f"| total {len(chosen)} | sha256 {digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
