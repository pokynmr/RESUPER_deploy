#!/usr/bin/python3
"""What the detected states actually differ in, for a multi-domain bundle.

Manuscript Section 3.9 uses this on 7SA5 (two-state apo Pin1): RESUPER finds the
WW and PPIase domains as two rigid ranges, the detector splits the 20 conformers
into two states, and this script asks whether the two states differ in the
*relative arrangement of the two domains*.  Both measures are computed per
conformer in its own frame, so they are invariant to overall rotation and
translation and comparable across conformers directly:

  separation  -- distance between the two domain Cα centroids;
  angle       -- angle between the two domains' principal axes (0..90 deg).

The per-state mean and spread are printed next to the whole-bundle figure. A
state difference that shows up here but not in the per-domain RMSD means the
states are different arrangements of the same rigid bodies -- which is exactly
what the residue-range method is meant to leave out of a single superposition.

    /usr/bin/python3 scripts/states_motion.py /tmp/7SA5.pdb

Read-only: prints a table, writes nothing.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from resuper.core import Params, identify_domains                      # noqa: E402
from resuper.pdb import parse_pdb                                      # noqa: E402
from resuper.states import detect_states                               # noqa: E402


def domain_ca(structure, conformer, chain, lo, hi):
    """Cα coordinates of one conformer over the reported range `lo..hi`."""
    return np.array([conformer[(r.chain, r.num)]["CA"]
                     for r in structure.residues
                     if (not chain or r.chain == chain) and lo <= r.num <= hi])


def measures(structure, ranges):
    """Per-conformer (separation, angle) for the first two reported ranges."""
    (c1, (l1, h1)), (c2, (l2, h2)) = ranges[:2]
    sep, ang = [], []
    for conf in structure.conformers:
        a = domain_ca(structure, conf, c1, l1, h1)
        b = domain_ca(structure, conf, c2, l2, h2)
        sep.append(float(np.linalg.norm(a.mean(axis=0) - b.mean(axis=0))))
        va = np.linalg.svd(a - a.mean(axis=0))[2][0]
        vb = np.linalg.svd(b - b.mean(axis=0))[2][0]
        ang.append(float(np.degrees(np.arccos(np.clip(abs(np.dot(va, vb)),
                                                      0.0, 1.0)))))
    return np.array(sep), np.array(ang)


def report(path):
    structure = parse_pdb(path)
    domains = identify_domains(structure, Params())
    if len(domains) < 2:
        print("%s: %d domain(s) -- interdomain arrangement needs at least two "
              "domains, skipped" % (os.path.basename(path), len(domains)))
        return
    pairs = []
    for d in domains:
        chains = d.chains or [""] * len(d.ranges)
        pairs += list(zip(chains, d.ranges))
    states = detect_states(structure)
    sep, ang = measures(structure, pairs)
    (c1, (l1, h1)), (c2, (l2, h2)) = pairs[:2]
    label = "%s%d..%s%d vs %s%d..%s%d" % (c1, l1, c1, h1, c2, l2, c2, h2)
    print("%s  (%d conformers, %d states)  domains %s"
          % (os.path.basename(path), structure.n_models, states.n_states,
             label))
    for name, values in (("separation (A)", sep), ("angle (deg)", ang)):
        print("  %-16s all: mean %6.2f sd %5.2f | %s"
              % (name, values.mean(), values.std(), "; ".join(
                  "state %d: mean %6.2f sd %5.2f"
                  % (i + 1, values[states.labels == i].mean(),
                     values[states.labels == i].std())
                  for i in range(states.n_states))))
    print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    for arg in sys.argv[1:]:
        report(arg)
