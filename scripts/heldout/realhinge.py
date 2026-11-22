"""Test set R: hinge motions applied to real NMR ensembles (real noise and disorder).

Pre-registered test set of scripts/heldout/PROTOCOL.md, committed after the candidate of
commit 580d0db was frozen and before it was run on the set.

Base ensembles: the first 20 (or all) models of twelve NMR ensembles that CYRANGE reports
as one domain and that are not used in any other held-out set: 1a9v, 1aca, 1adn, 1AYG,
1BA9, 1BM6, 1BZF, 1EH2, 1GM1, 1HZL, 1LKN, 1YZC (76-173 residues). Every model keeps its
own coordinates, so the internal disorder of the experimental ensemble is preserved; in
each model the motions of scripts/heldout/heldout.py are applied in torsion space (hinge:
phi and psi of one residue by N(0, A), A drawn from 10-20 deg; linker: phi and psi of
each residue by N(0, 30 deg)).

Scenarios per protein (positions as fractions of the length from a seeded generator):
  single   unchanged ensemble
  two      hinge at 0.3-0.7
  three    hinges at 0.25-0.4 and 0.6-0.75
  linker   disordered linker of 4-8 residues starting at 0.4-0.55
(48 ensembles) and the two-body scenario at A = 2, 4, 6, 8, 12, 16, 24 deg (84 ensembles).

Ground truth: hinge residues are not scored; linker residues are flexible, and so are
residues disordered in the experimental ensemble (CA RMSF > 2.0 A after iterative
superposition on the better-ordered half of the CA atoms); the remaining residues between hinges and linkers form
rigid bodies, those with at least 8 such residues the true domains. Scores as in
heldout.py.

usage (from the repo root):
  /usr/bin/python3 scripts/heldout/realhinge.py build     # -> cywork/realhinge/
  /usr/bin/python3 scripts/heldout/realhinge.py capture   # -> testdata/oracle_realhinge/
"""
import os, subprocess, sys
sys.path.insert(0, os.path.abspath("scripts/heldout"))
import numpy as np
import heldout as H

PROTEINS = ["1a9v", "1aca", "1adn", "1AYG", "1BA9", "1BM6", "1BZF", "1EH2", "1GM1",
            "1HZL", "1LKN", "1YZC"]
OUT = "cywork/realhinge"
REF = "testdata/oracle_realhinge"
REMOTE = "/tmp/resuper_realhinge"
NMAX = 20
FLEX_RMSF = 2.0
SWEEP = [2, 4, 6, 8, 12, 16, 24]


def read_models(code):
    """Per model: ATOM lines (first altloc only); all models share the atom list."""
    models, cur = [], None
    for line in open(f"testdata/{code}.pdb"):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = None
        elif cur is not None and line.startswith("ATOM") and line[16] in " A":
            cur.append(line)
    key = [l[12:27] for l in models[0]]
    models = [m for m in models if [l[12:27] for l in m] == key]
    return models[:NMAX]


def flexible(code):
    """Residue indices with CA RMSF > FLEX_RMSF in the experimental ensemble, after
    superposition on its better-ordered half (iterated: fit on the residues with RMSF
    below the median, recompute), so that disordered tails do not inflate the RMSF of
    the ordered core."""
    X = np.array([[[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                   for l in m if l[12:16] == " CA "] for m in read_models(code)])
    fit = np.arange(X.shape[1])
    for _ in range(10):
        ref = X[:, fit].mean(axis=0)
        for k in range(len(X)):
            P, Q = X[k, fit] - X[k, fit].mean(axis=0), ref - ref.mean(axis=0)
            U, _, Vt = np.linalg.svd(P.T @ Q)
            d = np.sign(np.linalg.det(U @ Vt))
            X[k] = (X[k] - X[k, fit].mean(axis=0)) @ (U @ np.diag([1, 1, d]) @ Vt) \
                + ref.mean(axis=0)
        rmsf = np.sqrt(((X - X.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))
        fit = np.where(rmsf <= np.median(rmsf))[0]
    return {i for i, v in enumerate(rmsf) if v > FLEX_RMSF}


def design():
    cases = {}
    for p, code in enumerate(PROTEINS):
        n = sum(1 for l in read_models(code)[0] if l[12:16] == " CA ")
        rng = np.random.default_rng(4000 + p)
        pos = lambda lo, hi: int(round(rng.uniform(lo, hi) * n))
        amp = lambda: float(rng.uniform(10, 20))
        h2 = (pos(0.3, 0.7), amp())
        h3 = [(pos(0.25, 0.4), amp()), (pos(0.6, 0.75), amp())]
        lk0, lkn = pos(0.4, 0.55), int(rng.integers(4, 9))
        scen = {"single": ([], []), "two": ([h2], []), "three": (h3, []),
                "linker": ([], list(range(lk0, lk0 + lkn)))}
        for s, (hinges, flex) in scen.items():
            cases[f"{code}_{s}"] = (code, hinges, flex, "exp")
        for a in SWEEP:
            cases[f"{code}_sweep{a:02d}"] = (code, [(h2[0], float(a))], [], "exp")
    return cases


def build_one(name, case, seed):
    code, hinges, flex, _ = case
    models = read_models(code)
    lines = models[0]
    names = [l[12:16].strip() for l in lines]
    keys = [(l[21], l[22:27]) for l in lines]
    order = list(dict.fromkeys(keys))
    ri = np.array([order.index(k) for k in keys])
    resn = {ri[i]: lines[i][17:20] for i in range(len(lines))}
    atom = {(ri[i], names[i]): i for i in range(len(lines))}
    moves = []
    for k, sd in hinges + [(k, H.TAIL_SD) for k in flex]:
        if resn[k] != "PRO":
            m = (ri > k) | ((ri == k) & np.array([n not in H.N_TERM_H for n in names]))
            moves.append((m, atom[(k, "N")], atom[(k, "CA")], sd))
        m = (ri > k) | ((ri == k) & np.array([n in H.O_ATOMS for n in names]))
        moves.append((m, atom[(k, "CA")], atom[(k, "C")], sd))
    rng = np.random.default_rng(seed)
    with open(f"{OUT}/{name}.pdb", "w") as f:
        f.write(f"TITLE     real-ensemble hinge test {name}\n")
        for c, model in enumerate(models):
            X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in model])
            for m, a, b, sd in moves:
                H.rotate(X, m, a, b, rng.normal(0, sd))
            f.write(f"MODEL     {c + 1:4d}\n")
            for l, p in zip(model, X):
                f.write(l[:16] + " " + l[17:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:])
            f.write("ENDMDL\n")
        f.write("END\n")


def truth(case):
    """heldout.truth with the experimentally disordered residues made flexible."""
    code, hinges, flex, _ = case
    disordered = flexible(code)
    lines, _ = H.read_base(code)
    nums = [int(l[22:26]) for l in lines if l[12:16] == " CA "]
    cut = sorted({k for k, _ in hinges} | set(flex))
    out, d, cur, runs = {}, 0, [], []
    for i in range(len(nums)):
        if i in cut:
            if cur:
                runs.append(cur)
            cur = []
        else:
            cur.append(i)
    if cur:
        runs.append(cur)
    for run in runs:
        rigid = [i for i in run if i not in disordered]
        if len(rigid) >= H.MU:
            for i in rigid:
                out[nums[i]] = d
            d += 1
        else:
            for i in rigid:
                out[nums[i]] = -1
        for i in run:
            if i in disordered:
                out[nums[i]] = -1
    for k in flex:
        out[nums[k]] = -1
    for k, _ in hinges:
        out[nums[k]] = None
    return out, d


def main():
    cases = design()
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "build":
        for s, (name, case) in enumerate(cases.items()):
            build_one(name, case, 11007 + s)
        print(f"{len(cases)} ensembles in {OUT}/")
    elif sys.argv[1] == "capture":
        os.makedirs(REF, exist_ok=True)
        subprocess.run(["ssh", H.ORACLE, f"rm -rf {REMOTE} && mkdir -p {REMOTE}"], check=True)
        subprocess.run(["scp", "-q"] + [f"{OUT}/{n}.pdb" for n in cases]
                       + [f"{H.ORACLE}:{REMOTE}/"], check=True)
        subprocess.run(["ssh", H.ORACLE, f"cd {REMOTE} && ls *.pdb | xargs -P 32 -I{{}} sh -c "
                        "'/usr/local/bin/cyrange {} > $(basename {} .pdb).out 2>&1; true'"],
                       check=True)
        tar = subprocess.run(["ssh", H.ORACLE, f"cd {REMOTE} && tar cf - *.out"],
                             capture_output=True, check=True).stdout
        subprocess.run(["tar", "xf", "-", "-C", REF], input=tar, check=True)
        print(f"{len(cases)} reference reports in {REF}/")


if __name__ == "__main__":
    main()
