"""Held-out ground-truth benchmark: eight proteins, torsion-space motions, realistic noise.

Pre-registered design. This file was committed before any method was run on the
ensembles it builds, and no resuper parameter may be changed in response to the results.
The sound-mode thresholds were set on ubiquitin (testdata/final.pdb) with rigid rotations
about a hinge atom and 0.05 A jitter (scripts/synthetic/truthbench.py); everything below
differs from that: the proteins, the motion model, the noise model and the design.

Base structures: model 1 of eight NMR ensembles, 110-222 residues, single chain, no gaps:
1A5J, 1a63, 1MZK, 1G03, 1INZ, 2GO9, 6R0J, 2I3E. 20 conformers per ensemble.

Motions are applied in torsion space, so bond lengths and angles are kept:
- a hinge at residue k perturbs phi_k and psi_k by N(0, A) in every conformer (A per hinge
  drawn uniformly from 10-20 deg; proline: psi only), moving everything C-terminal of it
  as one rigid body;
- a disordered tail or linker perturbs phi and psi of each of its residues by N(0, 30 deg).
Noise: every atom gets isotropic Gaussian noise of 0.05 A ("clean"). "Real" noise adds a
displacement of every residue as a whole, smoothed along the sequence (Gaussian kernel, SD 2
residues) so that neighbouring residues move together, with an SD of 0.3 x the CA RMSF of
the residue in the base NMR ensemble, clipped to 0.05-0.5 A.

Scenarios (chain positions as fractions of the length, drawn per protein from a seeded
generator):
  single   one body
  two      hinge at 0.3-0.7
  three    hinges at 0.25-0.4 and 0.6-0.75
  tail2    two bodies (hinge at 0.4-0.7) and a disordered N- or C-terminal tail of 8-16
  linker   two bodies joined by a disordered linker of 4-8 residues starting at 0.4-0.55
  tails1   one body with disordered tails of 8-16 residues at both ends
each with clean and real noise (96 ensembles), plus the two-body scenario with real noise
at A = 2, 4, 6, 8, 12, 16, 24 deg (56 ensembles; reported as a detection curve, since a
hinge of a few degrees has no right answer).

Ground truth: hinge residues are not scored; tail and linker residues are flexible; the
remaining runs of residues are rigid bodies, those of at least 8 residues the true domains
and shorter ones flexible.

Scores per method (resuper cyrange / sound, and the CYRANGE default report):
  count    number of predicted domains equals the number of true domains     (primary)
  mixing   a predicted domain draws at least 3 residues and 10% of its rigid residues
           from a second true domain                                          (primary)
  contam   fraction of the selected scored residues that are flexible         (primary)
  jaccard  mean over true domains of the best Jaccard with a predicted domain
  recall   fraction of true-domain residues selected

Set 2 (`--set 2`) repeats the design on eight further proteins with new seeds: 1a7m, 1DGQ,
1JBH, 2ARF, 1LM0, 1A67, 9Z30, 1D9S (101-188 residues). It was committed before any method
change made after set 1 was scored (the disordered-fragment test) was run on it.

Set 3 (`--set 3`) repeats the design again on 1JI8, 1I17, 1BR0, 2AXL, 1a5e, 1a5r, 1E8E,
1EII (103-156 residues); it is the first test set of scripts/heldout/PROTOCOL.md, committed
after the candidate of commit 580d0db was frozen and before it was run on the set.

usage (from the repo root):
  /usr/bin/python3 scripts/heldout/heldout.py build   [--set 2]  # ensembles -> cywork/heldout[2]/
  /usr/bin/python3 scripts/heldout/heldout.py capture [--set 2]  # reference -> testdata/oracle_heldout[2]/
  /usr/bin/python3 scripts/heldout/heldout.py score   [--set 2]
"""
import math, os, subprocess, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.abspath("scripts"))
import numpy as np
from eval_methods import parse_report
from resuper.cli import format_report
from resuper.core import Params, identify_domains
from resuper.pdb import parse_pdb

SETS = {
    "1": (["1A5J", "1a63", "1MZK", "1G03", "1INZ", "2GO9", "6R0J", "2I3E"], 1000, 7, ""),
    "2": (["1a7m", "1DGQ", "1JBH", "2ARF", "1LM0", "1A67", "9Z30", "1D9S"], 2000, 5007, "2"),
    "3": (["1JI8", "1I17", "1BR0", "2AXL", "1a5e", "1a5r", "1E8E", "1EII"], 3000, 9007, "3"),
}
ORACLE = os.environ.get("RESUPER_CYRANGE_HOST", "")  # ssh host with a licensed CYRANGE binary


def use_set(name):
    """Select a held-out set: its proteins, seeds and directories."""
    global SET, PROTEINS, DESIGN_SEED, BUILD_SEED, OUT, REF, REMOTE
    SET = name
    PROTEINS, DESIGN_SEED, BUILD_SEED, sfx = SETS[name]
    OUT = "cywork/heldout" + sfx
    REF = "testdata/oracle_heldout" + sfx
    REMOTE = "/tmp/resuper_heldout" + sfx


use_set(sys.argv[sys.argv.index("--set") + 1] if "--set" in sys.argv else "1")
NCONF = 20
TAIL_SD = 30.0
SWEEP = [2, 4, 6, 8, 12, 16, 24]
MU = 8
N_TERM_H = {"N", "H", "HN", "H1", "H2", "H3", "CA"}
O_ATOMS = {"O", "OXT", "OT1", "OT2"}


def read_base(code):
    """Model-1 ATOM lines and the per-residue CA RMSF of the whole ensemble."""
    lines, cas, cur = [], [], None
    for line in open(f"testdata/{code}.pdb"):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            cas.append(cur)
            cur = None
        elif line.startswith("ATOM") and line[16] in " A":
            if len(cas) == 0:
                lines.append(line)
            if line[12:16] == " CA " and cur is not None:
                cur.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    X = np.array(cas)
    for _ in range(5):
        mean = X.mean(axis=0)
        for k in range(len(X)):
            P, Q = X[k] - X[k].mean(axis=0), mean - mean.mean(axis=0)
            U, _, Vt = np.linalg.svd(P.T @ Q)
            d = np.sign(np.linalg.det(U @ Vt))
            X[k] = P @ (U @ np.diag([1, 1, d]) @ Vt) + mean.mean(axis=0)
    rmsf = np.sqrt(((X - X.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))
    return lines, rmsf


def design():
    """{name: (protein, hinges [(res index, A)], flexible residue indices, noise)}."""
    cases = {}
    for p, code in enumerate(PROTEINS):
        lines, _ = read_base(code)
        n = sum(1 for l in lines if l[12:16] == " CA ")
        rng = np.random.default_rng(DESIGN_SEED + p)
        pos = lambda lo, hi: int(round(rng.uniform(lo, hi) * n))
        amp = lambda: float(rng.uniform(10, 20))
        tail = lambda: int(rng.integers(8, 17))
        h2 = (pos(0.3, 0.7), amp())
        h3 = [(pos(0.25, 0.4), amp()), (pos(0.6, 0.75), amp())]
        ht = (pos(0.4, 0.7), amp())
        L = tail()
        t_n = bool(rng.integers(0, 2))
        tflex = list(range(L)) if t_n else list(range(n - L, n))
        lk0, lkn = pos(0.4, 0.55), int(rng.integers(4, 9))
        l1, l2 = tail(), tail()
        scen = {
            "single": ([], []),
            "two": ([h2], []),
            "three": (h3, []),
            "tail2": ([ht], tflex),
            "linker": ([], list(range(lk0, lk0 + lkn))),
            "tails1": ([], list(range(l1)) + list(range(n - l2, n))),
        }
        for s, (hinges, flex) in scen.items():
            for noise in ("clean", "real"):
                cases[f"{code}_{s}_{noise}"] = (code, hinges, flex, noise)
        for a in SWEEP:
            cases[f"{code}_sweep{a:02d}"] = (code, [(h2[0], float(a))], [], "real")
    return cases


def rotate(X, mask, a, b, deg):
    """Rotate X[mask] by deg about the axis a->b through b."""
    k = (X[b] - X[a]) / np.linalg.norm(X[b] - X[a])
    t = math.radians(deg)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + math.sin(t) * K + (1 - math.cos(t)) * K @ K
    X[mask] = (X[mask] - X[b]) @ R.T + X[b]


def build_one(name, case, seed):
    code, hinges, flex, noise = case
    lines, rmsf = read_base(code)
    names = [l[12:16].strip() for l in lines]
    keys = [(l[21], l[22:27]) for l in lines]
    order = list(dict.fromkeys(keys))
    ri = np.array([order.index(k) for k in keys])
    resn = {ri[i]: lines[i][17:20] for i in range(len(lines))}
    atom = {(ri[i], names[i]): i for i in range(len(lines))}
    base = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
    nres = len(order)
    sigma = np.clip(0.3 * rmsf, 0.05, 0.5)
    w = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    w /= np.sqrt((w ** 2).sum())
    moves = []          # (mask, axis atom a, axis atom b, sd)
    for k, sd in hinges + [(k, TAIL_SD) for k in flex]:
        if resn[k] != "PRO":
            m = (ri > k) | ((ri == k) & np.array([n not in N_TERM_H for n in names]))
            moves.append((m, atom[(k, "N")], atom[(k, "CA")], sd))
        m = (ri > k) | ((ri == k) & np.array([n in O_ATOMS for n in names]))
        moves.append((m, atom[(k, "CA")], atom[(k, "C")], sd))
    rng = np.random.default_rng(seed)
    with open(f"{OUT}/{name}.pdb", "w") as f:
        f.write(f"TITLE     held-out benchmark {name}\n")
        for c in range(NCONF):
            X = base.copy()
            for m, a, b, sd in moves:
                rotate(X, m, a, b, rng.normal(0, sd))
            X += rng.normal(0, 0.05, X.shape)
            if noise == "real":
                z = rng.normal(0, 1, (nres + 12, 3))
                field = np.stack([np.convolve(z[:, j], w, mode="valid") for j in range(3)], 1)
                X += (field * sigma[:, None])[ri]
            f.write(f"MODEL     {c + 1:4d}\n")
            for l, p in zip(lines, X):
                f.write(l[:16] + " " + l[17:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:])
            f.write("ENDMDL\n")
        f.write("END\n")


def truth(case):
    """Per residue number: domain id, -1 (flexible) or None (hinge, not scored)."""
    code, hinges, flex, _ = case
    lines, _ = read_base(code)
    nums = [int(l[22:26]) for l in lines if l[12:16] == " CA "]
    lab = [0] * len(nums)
    cut = sorted({k for k, _ in hinges} | set(flex))
    for k in flex:
        lab[k] = -1
    for k, _ in hinges:
        lab[k] = None
    d, runs = 0, []
    cur = []
    for i in range(len(nums)):
        if i in cut:
            if cur:
                runs.append(cur)
            cur = []
        else:
            cur.append(i)
    if cur:
        runs.append(cur)
    out = {}
    for run in runs:
        if len(run) >= MU:
            for i in run:
                out[nums[i]] = d
            d += 1
        else:
            for i in run:
                out[nums[i]] = -1
    for i in range(len(nums)):
        if i in cut:
            out[nums[i]] = lab[i]
    return out, d


def score(labels, ndom, pred):
    scored = {r for r, l in labels.items() if l is not None}
    flex = {r for r, l in labels.items() if l == -1}
    true = [{r for r, l in labels.items() if l == d} for d in range(ndom)]
    pred = [p & scored for p in pred]
    mixing = False
    for p in pred:
        counts = sorted((len(p & t) for t in true), reverse=True)
        rigid = sum(counts)
        if len(counts) > 1 and counts[1] >= 3 and counts[1] >= 0.1 * rigid:
            mixing = True
    allp = set().union(*pred) if pred else set()
    jac = [max((len(t & p) / len(t | p) for p in pred), default=0.0) for t in true]
    alltrue = set().union(*true) if true else set()
    return {"count": len(pred) == ndom, "mixing": mixing,
            "contam": len(allp & flex) / len(allp) if allp else 0.0,
            "jaccard": float(np.mean(jac)) if jac else 1.0,
            "recall": len(allp & alltrue) / len(alltrue) if alltrue else 1.0,
            "n": len(pred)}


def residues(text):
    return [{i for _, i in d[0]} for d in parse_report(text)]


def run_resuper(name):
    st = parse_pdb(f"{OUT}/{name}.pdb")
    return {m: format_report(f"{name}.pdb", st.n_models, identify_domains(st, Params(method=m)))
            for m in ("reference", "default")}


def main():
    cases = design()
    os.makedirs(OUT, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else "score"
    if cmd == "build":
        for s, (name, case) in enumerate(cases.items()):
            build_one(name, case, BUILD_SEED + s)
        print(f"{len(cases)} ensembles in {OUT}/")
    elif cmd == "capture":
        os.makedirs(REF, exist_ok=True)
        subprocess.run(["ssh", ORACLE, f"rm -rf {REMOTE} && mkdir -p {REMOTE}"], check=True)
        subprocess.run(["scp", "-q"] + [f"{OUT}/{n}.pdb" for n in cases]
                       + [f"{ORACLE}:{REMOTE}/"], check=True)
        subprocess.run(["ssh", ORACLE, f"cd {REMOTE} && ls *.pdb | xargs -P 32 -I{{}} sh -c "
                        "'/usr/local/bin/cyrange {} > $(basename {} .pdb).out 2>&1; true'"], check=True)
        tar = subprocess.run(["ssh", ORACLE, f"cd {REMOTE} && tar cf - *.out"],
                             capture_output=True, check=True).stdout
        subprocess.run(["tar", "xf", "-", "-C", REF], input=tar, check=True)
        print(f"{len(cases)} reference reports in {REF}/")
    else:
        names = list(cases)
        with Pool() as pool:
            ours = dict(zip(names, pool.map(run_resuper, names)))
        methods = ["baseline", "reference", "default"]
        rows = {}
        for name in names:
            labels, ndom = truth(cases[name])
            rep = {"baseline": open(f"{REF}/{name}.out").read(), **ours[name]}
            rows[name] = (ndom, {m: score(labels, ndom, residues(rep[m])) for m in methods})
        main_names = [n for n in names if "sweep" not in n]
        groups = {"all 96": main_names}
        for key in ("clean", "real"):
            groups[f"{key} noise"] = [n for n in main_names if n.endswith(key)]
        for s in ("single", "two", "three", "tail2", "linker", "tails1"):
            groups[s] = [n for n in main_names if f"_{s}_" in n]
        print(f"{'group':12s} {'n':>3s}  " + "  ".join(
            f"{m:>40s}" for m in methods))
        print(" " * 18 + "  ".join(f"{'count mix  contam  jacc  recall':>40s}" for _ in methods))
        for g, ns in groups.items():
            cells = []
            for m in methods:
                r = [rows[n][1][m] for n in ns]
                cells.append(f"{sum(x['count'] for x in r):5d} {sum(x['mixing'] for x in r):4d}  "
                             f"{np.mean([x['contam'] for x in r]):6.3f} "
                             f"{np.mean([x['jaccard'] for x in r]):5.3f}  "
                             f"{np.mean([x['recall'] for x in r]):6.3f}")
            print(f"{g:12s} {len(ns):3d}  " + "  ".join(f"{c:>40s}" for c in cells))
        print("\ntwo-body detection (ensembles with 2 predicted domains, out of 8) by hinge SD:")
        print("deg  " + "  ".join(f"{m:>9s}" for m in methods))
        for a in SWEEP:
            ns = [n for n in names if n.endswith(f"sweep{a:02d}")]
            print(f"{a:3d}  " + "  ".join(
                f"{sum(rows[n][1][m]['n'] == 2 for n in ns):9d}" for m in methods))
        print("\nper-case misses (domain count wrong or mixing):")
        for name in main_names:
            ndom, r = rows[name]
            bad = [f"{m}:n{r[m]['n']}{'M' if r[m]['mixing'] else ''}" for m in methods
                   if not r[m]["count"] or r[m]["mixing"]]
            if bad:
                print(f"  {name:22s} true {ndom}  " + "  ".join(bad))


if __name__ == "__main__":
    main()
