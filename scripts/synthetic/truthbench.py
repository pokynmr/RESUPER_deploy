"""Ground-truth domain benchmark on synthetic intact-hinge ensembles.

Each case is a segments.py spec on testdata/final.pdb (76 residues, 10 models, every
residue jittered 0.05 A). Residues that receive exactly the same set of applied
rotations move as one rigid body; bodies of >= 8 residues are the true domains, smaller
independently moving pieces (bending tails, flexible linkers) belong to no domain, and
the two residues at each pivot are ambiguous and not scored.

Scores per method (resuper --method reference / sound, and the reference's default report
when captured): domain count correct; mean Jaccard of every true domain with its best
matching predicted domain (unscored residues removed); contamination = fraction of the
predicted residues that are flexible (no-domain) residues.

usage (from the repo root):
  /usr/bin/python3 scripts/synthetic/truthbench.py --capture   # build + run the reference
  /usr/bin/python3 scripts/synthetic/truthbench.py             # score
Ensembles and reference reports go to cywork/truth/ (specs in TRUTH below).
"""
import os, re, sys
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segments
from resuper.core import Params, identify_domains
from resuper.pdb import parse_pdb

OUT = "cywork/truth"
J = 0.05


def body(lo, hi, deg, pivot):
    return f"{lo}-{hi}:{deg}:{J}:p{pivot}"


def tail_n(L, deg):
    """N-terminal tail 1..L bending at a hinge every 4 residues (nested pivots)."""
    return ",".join(f"1-{h - 1}:{deg}:{J}:p{h}" for h in range(L + 1, 4, -4))


def linker(lo, hi, deg):
    """Flexible linker lo..hi: every residue a hinge, nested toward the C-terminus."""
    return ",".join(f"{p}-76:{deg}:{J}:p{p}" for p in range(lo, hi + 1))


TRUTH = {
    "t1_single": f"1-76:0:{J}",
    "t2_h20_15": f"1-76:0:{J}," + body(20, 76, 15, 20),
    "t2_h38_15": f"1-76:0:{J}," + body(38, 76, 15, 38),
    "t2_h38_30": f"1-76:0:{J}," + body(38, 76, 30, 38),
    "t2_h55_20": f"1-76:0:{J}," + body(55, 76, 20, 55),
    "t3_eq_20": f"1-76:0:{J}," + body(26, 76, 20, 26) + "," + body(51, 76, 20, 51),
    "t3_nest_30": f"1-76:0:{J}," + body(26, 76, 20, 26) + "," + body(51, 76, 30, 51),
    "t3_mid12": f"1-76:0:{J}," + body(33, 76, 25, 33) + "," + body(45, 76, 25, 45),
    "tail12": f"1-76:0:{J}," + tail_n(12, 20),
    "tail16": f"1-76:0:{J}," + tail_n(16, 25),
    "tail16_2dom": f"1-76:0:{J}," + body(45, 76, 25, 45) + "," + tail_n(16, 25),
    "link6_20": f"1-76:0:{J}," + linker(36, 41, 20),
    "link4_30": f"1-76:0:{J}," + linker(37, 40, 30),
}


def truth_labels(spec, nres=76, mu=8):
    """Per residue: domain id (0..), or -1 (flexible), or None (pivot, unscored)."""
    segs, hinge = segments.parse_spec(spec)
    sig = {}
    for r in range(1, nres + 1):
        sig[r] = tuple(k for k, (lo, hi, deg, _) in enumerate(segs) if lo <= r <= hi and deg > 0)
    groups = {}
    for r, s in sig.items():
        groups.setdefault(s, []).append(r)
    labels, d = {}, 0
    for s, members in sorted(groups.items(), key=lambda kv: min(kv[1])):
        if len(members) >= mu:
            for r in members:
                labels[r] = d
            d += 1
        else:
            for r in members:
                labels[r] = -1
    for h in hinge:
        if h is not None:
            for r in (h - 1, h):
                if r in labels:
                    labels[r] = None
    return labels, d


def residue_sets(text):
    out = []
    for m in re.finditer(r"Optimal range (.*?): RMSD", text):
        s = set()
        for seg in m.group(1).split(","):
            a, _, b = seg.strip().partition("..")
            s |= set(range(int(a), int(b or a) + 1))
        out.append(s)
    return out


def score(labels, ndom, pred):
    scored = {r for r, l in labels.items() if l is not None}
    flex = {r for r, l in labels.items() if l == -1}
    true = [{r for r, l in labels.items() if l == d} for d in range(ndom)]
    pred = [p & scored for p in pred]
    jac = []
    for t in true:
        jac.append(max((len(t & p) / len(t | p) for p in pred), default=0.0))
    allp = set().union(*pred) if pred else set()
    contam = len(allp & flex) / len(allp) if allp else 0.0
    return len(pred) == ndom, sum(jac) / len(jac) if jac else 1.0, contam


def main():
    os.makedirs(OUT, exist_ok=True)
    if "--capture" in sys.argv:
        for name, spec in TRUTH.items():
            segments.generate(f"{OUT}/{name}.pdb", spec)
        segments.oracle([(n, s) for n, s in TRUTH.items()], OUT)
        return
    methods = ["baseline", "reference", "default"]
    tot = {m: [0, 0.0, 0.0, 0] for m in methods}
    print(f"{'case':14s} ntrue  " + "  ".join(f"{m:>22s}" for m in methods))
    for name, spec in TRUTH.items():
        labels, ndom = truth_labels(spec)
        row = []
        for meth in methods:
            if meth == "baseline":
                p = f"{OUT}/{name}.default"
                if not os.path.exists(p):
                    row.append(f"{'-':>22s}")
                    continue
                pred = residue_sets(open(p).read())
            else:
                st = parse_pdb(f"{OUT}/{name}.pdb")
                pred = [{i for lo, hi in d.ranges for i in range(lo, hi + 1)}
                        for d in identify_domains(st, Params(method=meth))]
            ok, jac, con = score(labels, ndom, pred)
            t = tot[meth]
            t[0] += ok; t[1] += jac; t[2] += con; t[3] += 1
            row.append(f"n{len(pred)}{'+' if ok else '-'} J{jac:.2f} c{con:.2f}".rjust(22))
        print(f"{name:14s} {ndom:5d}  " + "  ".join(row))
    print()
    for m, (ok, jac, con, n) in tot.items():
        if n:
            print(f"{m:10s} count correct {ok}/{n}  mean Jaccard {jac / n:.3f}  "
                  f"mean contamination {con / n:.3f}")


if __name__ == "__main__":
    main()
