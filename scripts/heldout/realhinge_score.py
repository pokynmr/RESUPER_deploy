"""Score the real-hinge set (12 experimental NMR ensembles, disorder preserved)
against the reference baseline. Uses the same truth/score as heldout.py.

usage: /usr/bin/python3 scripts/heldout/realhinge_score.py
"""
import os, sys
sys.path.insert(0, os.path.abspath("scripts/heldout"))
sys.path.insert(0, os.path.abspath("scripts"))
import realhinge as RH
from eval_methods import parse_report
from heldout import truth, score
from resuper.cli import format_report
from resuper.core import Params, identify_domains
from resuper.pdb import parse_pdb

OUT = "cywork/realhinge"
REF = "testdata/oracle_realhinge"


def residues(text):
    return [{i for _, i in d[0]} for d in parse_report(text)]


def all_names():
    return list(RH.design().keys())


def run_all(n, st):
    """Format every method for one ensemble; parse the PDB exactly once."""
    ndom, labels = truth(RH.design()[n])
    rep = {m: format_report(f"{n}.pdb", st.n_models,
                            identify_domains(st, Params(method=m)))
           for m in ("reference", "default")}
    rep["baseline"] = open(f"{REF}/{n}.out").read()
    return ndom, labels, rep


def main():
    alln = all_names()
    groups = {"all": alln}
    for s in ("single", "two", "three", "linker"):
        groups[s] = [n for n in alln if f"_{s}" in n]
    print(f"{'group':12s} {'n':>3s}  " + "  ".join(
        f"{m:>40s}" for m in ("baseline", "reference", "default")))
    print(" " * 18 + "  ".join(f"{'count mix  contam  jacc  recall':>40s}" for _ in range(3)))
    for g, ns in groups.items():
        cells = []
        for m in ("baseline", "reference", "default"):
            r = []
            for n in ns:
                st = parse_pdb(f"{OUT}/{n}.pdb")
                ndom, labels, rep = run_all(n, st)
                r.append(score(labels, ndom, residues(rep[m])))
            cells.append(f"{sum(x['count'] for x in r):5d} {sum(x['mixing'] for x in r):4d}  "
                         f"{sum(x['contam'] for x in r)/len(r):6.3f} "
                         f"{sum(x['jaccard'] for x in r)/len(r):5.3f}  "
                         f"{sum(x['recall'] for x in r)/len(r):6.3f}")
        print(f"{g:12s} {len(ns):3d}  " + "  ".join(f"{c:>40s}" for c in cells))


if __name__ == "__main__":
    main()
