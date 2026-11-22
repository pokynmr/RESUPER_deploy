"""Input-perturbation analysis on the real evaluation ensembles.

Pre-registered design, committed before any method was run on the perturbed inputs.

A good range-selection method gives nearly the same answer when an ensemble is perturbed
in ways that do not change its structure. For each of the 134 real ensembles of
scripts/eval_methods.py (reference report with at least one range), six inputs are written:
  ident     the ensemble rewritten unchanged (baseline; must reproduce the original report)
  reverse   the models in reverse order (every method should be invariant)
  drop1     a random 80% of the models (seed 1)
  drop2     a random 80% of the models (seed 2)
  noise10   isotropic Gaussian noise of 0.10 A on every atom (seed 1)
  noise25   isotropic Gaussian noise of 0.25 A on every atom (seed 1)
Every method (CYRANGE, resuper cyrange, resuper sound) is compared with its own output on
"ident": the fraction of ensembles with the same number of domains, the mean Jaccard index
of the selected residue sets, and the fraction of identical range lists.

usage (from the repo root):
  /usr/bin/python3 scripts/heldout/perturb.py build     # inputs -> cywork/perturb/
  /usr/bin/python3 scripts/heldout/perturb.py capture   # reference -> testdata/oracle_perturb/
  /usr/bin/python3 scripts/heldout/perturb.py score
"""
import os, subprocess, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.abspath("scripts"))
import numpy as np
from eval_cases import CURATED, AISAR, pacsy_ids
from eval_methods import parse_report
from resuper.cli import format_report
from resuper.core import Params, identify_domains
from resuper.pdb import parse_pdb

OUT = "cywork/perturb"
REF = "testdata/oracle_perturb"
ORACLE, REMOTE = os.environ.get("RESUPER_CYRANGE_HOST", ""), "/tmp/resuper_perturb"
PERTS = ["ident", "reverse", "drop1", "drop2", "noise10", "noise25"]
METHODS = ["baseline", "reference", "default"]


def cases():
    out = []
    for base in dict.fromkeys(CURATED + AISAR + pacsy_ids()):
        o, p = f"testdata/{base}.out", f"testdata/{base}.pdb"
        if os.path.exists(o) and os.path.exists(p) and "Optimal range" in open(o).read():
            out.append(base)
    return out


def split_models(path):
    """(lines before the first model, [model body lines], lines after the last model)."""
    head, models, tail, cur = [], [], [], None
    for line in open(path):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = None
        elif cur is not None:
            cur.append(line)
        elif models:
            tail.append(line)
        else:
            head.append(line)
    return head, models, tail


def write(path, head, models, tail):
    with open(path, "w") as f:
        f.writelines(head)
        for i, body in enumerate(models):
            f.write(f"MODEL     {i + 1:4d}\n")
            f.writelines(body)
            f.write("ENDMDL\n")
        f.writelines(tail)


def noisy(body, rng, sd):
    out = []
    for l in body:
        if l.startswith(("ATOM", "HETATM")):
            x = [float(l[30:38]), float(l[38:46]), float(l[46:54])] + rng.normal(0, sd, 3)
            l = l[:30] + f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}" + l[54:]
        out.append(l)
    return out


def build(base):
    head, models, tail = split_models(f"testdata/{base}.pdb")
    n = len(models)
    keep = lambda seed: sorted(np.random.default_rng(seed).choice(n, round(0.8 * n),
                                                                   replace=False))
    variants = {
        "ident": models,
        "reverse": models[::-1],
        "drop1": [models[i] for i in keep(1)],
        "drop2": [models[i] for i in keep(2)],
        "noise10": [noisy(b, np.random.default_rng(1 + i), 0.10) for i, b in enumerate(models)],
        "noise25": [noisy(b, np.random.default_rng(1 + i), 0.25) for i, b in enumerate(models)],
    }
    for k, v in variants.items():
        write(f"{OUT}/{base}__{k}.pdb", head, v, tail)


def run_resuper(name):
    st = parse_pdb(f"{OUT}/{name}.pdb")
    return {m: format_report(f"{name}.pdb", st.n_models, identify_domains(st, Params(method=m)))
            for m in ("reference", "default")}


def ranges(text):
    return [d[0] for d in parse_report(text)]


def main():
    bases = cases()
    names = [f"{b}__{p}" for b in bases for p in PERTS]
    os.makedirs(OUT, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else "score"
    if cmd == "build":
        for b in bases:
            build(b)
        print(f"{len(names)} inputs from {len(bases)} ensembles in {OUT}/")
    elif cmd == "capture":
        os.makedirs(REF, exist_ok=True)
        subprocess.run(["ssh", ORACLE, f"rm -rf {REMOTE} && mkdir -p {REMOTE}"], check=True)
        subprocess.run(["scp", "-q"] + [f"{OUT}/{n}.pdb" for n in names]
                       + [f"{ORACLE}:{REMOTE}/"], check=True)
        subprocess.run(["ssh", ORACLE, f"cd {REMOTE} && ls *.pdb | xargs -P 32 -I{{}} sh -c "
                        "'/usr/local/bin/cyrange {} > $(basename {} .pdb).out 2>&1; true'"], check=True)
        tar = subprocess.run(["ssh", ORACLE, f"cd {REMOTE} && tar cf - *.out"],
                             capture_output=True, check=True).stdout
        subprocess.run(["tar", "xf", "-", "-C", REF], input=tar, check=True)
        print(f"{len(names)} reference reports in {REF}/")
    else:
        with Pool() as pool:
            ours = dict(zip(names, pool.map(run_resuper, names)))
        rep = {n: {"baseline": open(f"{REF}/{n}.out").read(), **ours[n]} for n in names}
        orig_ok = {m: 0 for m in METHODS}
        for b in bases:
            want = open(f"testdata/{b}.out").read()
            r = rep[f"{b}__ident"]
            orig_ok["baseline"] += ranges(r["baseline"]) == ranges(want)
            st = parse_pdb(f"testdata/{b}.pdb")
            for m in ("reference", "default"):
                got = format_report("x", st.n_models, identify_domains(st, Params(method=m)))
                orig_ok[m] += ranges(r[m]) == ranges(got)
        print(f"{len(bases)} ensembles; 'ident' reproduces the original ranges: "
              + ", ".join(f"{m} {orig_ok[m]}" for m in METHODS))
        print(f"\n{'perturbation':12s}" + "".join(f"{m:>30s}" for m in METHODS))
        print(" " * 12 + "".join(f"{'same n  Jaccard  identical':>30s}" for _ in METHODS))
        unstable = {m: {} for m in METHODS}
        for p in PERTS[1:]:
            cells = []
            for m in METHODS:
                same, jac, ident = 0, [], 0
                for b in bases:
                    a, c = ranges(rep[f"{b}__ident"][m]), ranges(rep[f"{b}__{p}"][m])
                    sa = set().union(*a) if a else set()
                    sc = set().union(*c) if c else set()
                    same += len(a) == len(c)
                    jac.append(len(sa & sc) / max(1, len(sa | sc)))
                    ident += a == c
                    if len(a) != len(c):
                        unstable[m].setdefault(b, []).append(f"{p}:{len(a)}>{len(c)}")
                cells.append(f"{same:6d}  {np.mean(jac):7.3f}  {ident:9d}")
            print(f"{p:12s}" + "".join(f"{c:>30s}" for c in cells))
        for m in METHODS:
            print(f"\n{m}: ensembles whose domain count changes under some perturbation "
                  f"({len(unstable[m])}):")
            for b, v in unstable[m].items():
                print(f"  {b:10s} {' '.join(v)}")


if __name__ == "__main__":
    main()
