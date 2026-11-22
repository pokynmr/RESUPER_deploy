"""Compare resuper methods against the reference reports on the evaluation set.

For every case of scripts/eval_cases.py with a usable oracle report, run each method
in-process and report, per method:
  exact    byte-identical report text (meaningful for --method reference only)
  ndom     same number of domains as the reference
  jaccard  mean Jaccard index of the selected residue sets (union over domains)
  cover    mean fraction of residues selected (reference in brackets)
  rmsd     mean residue-weighted RMSD over domains (reference in brackets)
and list the cases whose domain count differs from the reference.

usage: /usr/bin/python3 scripts/eval_methods.py [method ...]   (default: reference default)
"""
import os, re, sys
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from eval_cases import CURATED, AISAR, pacsy_ids
from resuper.cli import format_report
from resuper.core import Params, identify_domains
from resuper.pdb import parse_pdb

RANGE = re.compile(r"Optimal range (.*?): RMSD ([0-9.]+) A, \d+ gaps?, (\d+) residues?\.")


def residue_set(spec):
    out = set()
    for seg in spec.split(","):
        seg = seg.strip()
        m = re.match(r"([A-Za-z]?)(-?\d+)(?:\.\.[A-Za-z]?(-?\d+))?$", seg)
        if m:
            a, b = int(m.group(2)), int(m.group(3) or m.group(2))
            out |= {(m.group(1), i) for i in range(a, b + 1)}
    return out


def parse_report(text):
    doms = [(residue_set(m.group(1)), float(m.group(2)), int(m.group(3)))
            for m in RANGE.finditer(text)]
    return doms


def summary(doms):
    sel = set().union(*[d[0] for d in doms]) if doms else set()
    n = sum(d[2] for d in doms)
    rmsd = sum(d[1] * d[2] for d in doms) / n if n else float("nan")
    return sel, rmsd


def main():
    methods = sys.argv[1:] or ["reference", "default"]
    cases = list(dict.fromkeys(CURATED + AISAR + pacsy_ids()))
    stats = {m: {"exact": 0, "ndom": 0, "jac": [], "cov": [], "rmsd": [], "diff": []}
             for m in methods}
    ref_cov, ref_rmsd, n = [], [], 0
    for base in cases:
        out, pdb = f"testdata/{base}.out", f"testdata/{base}.pdb"
        if not (os.path.exists(out) and os.path.exists(pdb)):
            continue
        text = open(out).read()
        if "Optimal range" not in text:
            continue
        ref = parse_report(text)
        st = parse_pdb(pdb)
        nres = len(st.residues)
        rsel, rr = summary(ref)
        ref_cov.append(len(rsel) / nres)
        ref_rmsd.append(rr)
        n += 1
        for meth in methods:
            doms = identify_domains(st, Params(method=meth))
            rep = format_report(os.path.basename(pdb), st.n_models, doms)
            s = stats[meth]
            s["exact"] += rep.rstrip("\n") == text.rstrip("\n")
            got = parse_report(rep)
            s["ndom"] += len(got) == len(ref)
            if len(got) != len(ref):
                s["diff"].append(f"{base}({len(got)}/{len(ref)})")
            gsel, gr = summary(got)
            s["jac"].append(len(gsel & rsel) / max(1, len(gsel | rsel)))
            s["cov"].append(len(gsel) / nres)
            s["rmsd"].append(gr)
    print(f"{n} cases; reference: cover {np.mean(ref_cov):.3f}, rmsd {np.nanmean(ref_rmsd):.3f}")
    for meth in methods:
        s = stats[meth]
        print(f"{meth:8s} exact {s['exact']:3d}/{n}  ndom {s['ndom']:3d}/{n}  "
              f"jaccard {np.mean(s['jac']):.3f}  cover {np.mean(s['cov']):.3f}  "
              f"rmsd {np.nanmean(s['rmsd']):.3f}")
        print(f"   domain-count differences: {' '.join(s['diff'])}")


if __name__ == "__main__":
    main()
