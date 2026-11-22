#!/usr/bin/python3
"""RMSD of RESUPER's reported ranges, over the whole bundle and per detected state.

Manuscript Section 3.9 uses this to ask what a detected split is worth in the
quantity the method actually reports. It is not a test of whether the states are
real: when a bundle's states are two arrangements of the same rigid bodies, each
range is as variable inside a state as across the bundle, and the RMSD hardly
moves (7SA5 0.92 -> 0.94 and 0.87). Where a second state is not just a
re-arrangement it does respond (8SG2 0.58 -> 0.37 and 0.60). That contrast is why
the state count is reported separately instead of being inferred from the RMSD.

The RMSD is exactly the quantity RESUPER reports -- rmsd_to_mean_converged over
the N/CA/C' atoms of the range, i.e. the mean per-conformer distance to the mean
after optimal mutual superposition -- computed here on a subset of the
conformers instead of all of them.

    PYTHONPATH=src /usr/bin/python3 scripts/states_rmsd.py \
        testdata/7sa5.pdb testdata/7suq.pdb testdata/8sg2.pdb

Read-only: prints a table per structure, writes nothing.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from resuper.core import Params, _backbone_matrix, identify_domains   # noqa: E402
from resuper.geometry import rmsd_to_mean_converged                   # noqa: E402
from resuper.pdb import parse_pdb                                     # noqa: E402
from resuper.states import detect_states                              # noqa: E402


def range_residues(structure, chain, lo, hi):
    """0-based residue indices of the reported range `lo..hi` on `chain`."""
    return [i for i, r in enumerate(structure.residues)
            if (not chain or r.chain == chain) and lo <= r.num <= hi]


def range_rmsd(structure, conf_indices, chain, lo, hi):
    """RESUPER's reported RMSD for one range, restricted to given conformers."""
    rng = range_residues(structure, chain, lo, hi)
    if not rng:
        return None, 0
    ens, _ = _backbone_matrix(structure, [structure.conformers[i]
                                          for i in conf_indices], rng,
                              exclude_cterm=False)
    return rmsd_to_mean_converged(ens), len(rng)


def report(path):
    structure = parse_pdb(path)
    domains = identify_domains(structure, Params())
    states = detect_states(structure)
    all_idx = list(range(structure.n_models))
    print("%s  (%d conformers, %d state%s)"
          % (os.path.basename(path), structure.n_models, states.n_states,
             "" if states.n_states == 1 else "s"))
    header = "  %-14s %8s %8s" % ("range", "residues", "all")
    for s in range(states.n_states):
        header += " %10s" % ("state %d" % (s + 1))
    print(header)
    for d in domains:
        chains = d.chains or [""] * len(d.ranges)
        for chain, (lo, hi) in zip(chains, d.ranges):
            label = "%s%d..%s%d" % (chain, lo, chain, hi) if chain \
                else "%d..%d" % (lo, hi)
            rms_all, n = range_rmsd(structure, all_idx, chain, lo, hi)
            row = "  %-14s %8d %8.2f" % (label, n, rms_all)
            for st in states.states:
                rms, _ = range_rmsd(structure, st.indices, chain, lo, hi)
                row += " %10.2f" % rms
            print(row)
    print("  state membership: " + "; ".join(
        "state %d = %s" % (i + 1, ", ".join(str(m) for m in st.members))
        for i, st in enumerate(states.states)))
    print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    for arg in sys.argv[1:]:
        report(arg)
