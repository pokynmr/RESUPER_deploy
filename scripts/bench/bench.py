"""Benchmark B: domain detection, disorder, multi-state ensembles and input errors.

Compares the CYRANGE binary (baseline), resuper --method reference and resuper
--method default on 40 solution-NMR proteins sampled from the PDB (scripts/bench/pool.py,
scripts/bench/proteins.txt). The design was committed before any method was run on it.

Motions are applied in torsion space as in scripts/heldout/heldout.py (rotation about
the N-CA bond for phi and the CA-C bond for psi; proline: psi only). "Synthetic" bases
clone model 1 and add the realistic noise of heldout.py (whole-residue displacements
smoothed along the sequence, SD 0.3 x CA RMSF of the residue in the NMR ensemble after
superposition on its better-ordered half, clipped to 0.05-0.5 A, plus 0.05 A per atom);
"experimental" bases use the deposited models. 20 conformers unless stated.

Suites (n = chain length; positions drawn per protein from a seeded generator):
  clusters  synthetic base, k rigid bodies joined by hinges (phi, psi ~ N(0, A),
            A 10-20 deg): c1 (none), c2 (hinge at 0.3-0.7 n), c3 (0.2-0.4, 0.6-0.8),
            c4 (0.15-0.3, 0.4-0.6, 0.7-0.85)
  expclus   experimental base, only proteins the reference reports as one domain:
            e1, e2, e3 as c1-c3; residues with CA RMSF > 2 A in the deposited ensemble
            are flexible
  idr       synthetic base: idr_term (N- or C-terminal 25-40% of n disordered),
            idr_both (both termini, 15-25% each), idr_loop (internal disordered
            region of 15-25 residues starting at 0.35-0.5 n); disorder = phi, psi
            ~ N(0, 40 deg) and every side-chain dihedral ~ N(0, 60 deg) per residue
            and conformer (v2; v1 left the side chains rigid, so disordered
            residues kept order parameters of 1 and counted as core atoms)
  idp       synthetic base, every residue disordered (idp_full; no true domain),
            except one 12-residue segment at 0.3-0.6 n (idp_seg; one domain) or two
            12-residue segments at 0.15-0.3 n and 0.6-0.75 n (idp_2seg; two domains)
  multi     synthetic base, discrete states of a hinge at 0.35-0.65 n (phi and psi of
            the hinge residue changed by +theta in state B, -theta in state C, plus
            N(0, 3 deg) per conformer): ms2_15 and ms2_30 (states A/B 10:10, theta 15
            or 30 deg), ms2_skew (16:4, 30 deg), ms3 (A/B/C 7:7:6, 25 deg); and a
            two-body hinge (A 10-20 deg) with nc5, nc10, nc40 conformers. Two true
            domains throughout.
  errors    defective inputs made from the deposited ensembles of the first 8 pool
            proteins (DEFECTS below), each compared with the same method's output on
            the unchanged file ("clean").

Ground truth as in heldout.py: hinge residues are not scored; disordered residues are
flexible; runs of rigid residues between cuts with at least 8 residues are true
domains, shorter runs flexible.

Scores per method and scenario: domain count correct; mixed domains; adjusted Rand
index of the residue partition (true domain or flexible against predicted domain or
unselected, over scored residues); contamination (flexible residues among the selected);
Jaccard of each true domain with its best predicted domain; for idp_full the number of
reported domains and the fraction of residues selected. Errors: outcome per input,
crash (traceback, segmentation fault, timeout), refused (non-zero exit with a message),
near (same domain count as clean and Jaccard >= 0.8), drift (output that differs) or
empty (exit 0, no range); for defects without a sensible answer any non-crash outcome
is acceptable.

usage (from the repo root; add `--bench C` for benchmark C):
  /usr/bin/python3 scripts/bench/bench.py build-base # deposited ensembles -> cywork/bench/in/
  /usr/bin/python3 scripts/bench/bench.py capture base_<id> ...   # (expclus needs these)
  /usr/bin/python3 scripts/bench/bench.py build      # -> cywork/bench/in/
  /usr/bin/python3 scripts/bench/bench.py capture    # reference -> testdata/oracle_bench/
  /usr/bin/python3 scripts/bench/bench.py run [prefix ...]   # resuper -> cywork/bench/out/
  /usr/bin/python3 scripts/bench/bench.py run [prefix ...] KEY=VALUE ... tag=<t>
        # default method with resuper.core constants overridden -> cywork/bench/out/default<t>/
  /usr/bin/python3 scripts/bench/bench.py score [tag=<t>] [suites=clusters,idp,errors,...]
"""
import ast, json, math, os, subprocess, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.abspath("scripts"))
sys.path.insert(0, os.path.abspath("scripts/heldout"))
import numpy as np
import heldout as H
from eval_methods import parse_report

PDB = "cywork/bench/pdb"
# Benchmark B (default) or C (`--bench C`: the test set of round 2 of PROTOCOL.md, same
# design on the proteins of proteins_c.txt with other seeds).
BSET = sys.argv[sys.argv.index("--bench") + 1] if "--bench" in sys.argv else "B"
_SFX = "" if BSET == "B" else BSET
IN = f"cywork/bench{_SFX}/in"
OUT = f"cywork/bench{_SFX}/out"
REF = f"testdata/oracle_bench{_SFX}"
REMOTE = f"/tmp/resuper_bench{_SFX}"
LIST = "scripts/bench/proteins.txt" if BSET == "B" else f"scripts/bench/proteins_{BSET.lower()}.txt"
SEED_OFFSET = 0 if BSET == "B" else 1000
MU = 8
FLEX_RMSF = 2.0
IDR_SD = 40.0
CHI_SD = 60.0
TIMEOUT = 120
METHODS = ["baseline", "reference", "default"]
TAG = ""        # score: read the sound outputs of a configuration run (run ... tag=<t>)
SUITES = {"clusters": ["c1", "c2", "c3", "c4"], "expclus": ["e1", "e2", "e3"],
          "idr": ["idr_term", "idr_both", "idr_loop"],
          "idp": ["idp_full", "idp_seg", "idp_2seg"],
          "multi": ["ms2_15", "ms2_30", "ms2_skew", "ms3", "nc5", "nc10", "nc40"]}
# defect: expectation ("near": a near-clean answer exists; "any": no sensible answer)
DEFECTS = {
    "clean": "near", "miss_sidechain_atom": "near", "miss_ca_all": "near",
    "miss_ca_one_model": "near", "gap3": "near", "hetatm_ligand": "near", "mse": "near",
    "altloc": "near", "icode": "near", "strip_h": "near", "split_chain": "near",
    "dup_atom": "near", "miss_residue_one_model": "near", "shuffled_atoms": "near",
    "crlf": "near", "model_numbers": "near", "far_model": "near", "two_models": "any",
    "single_model": "any", "bare_single_model": "any", "identical_models": "any",
    "nan_coord": "any", "junk_lines": "any", "empty": "any", "no_model_records": "any",
    "tiny": "any",
}


def proteins():
    return open(LIST).read().split()


# ---------------------------------------------------------------- base structures

def read_models(code):
    models, cur = [], None
    for line in open(f"{PDB}/{code}.pdb"):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = None
        elif cur is not None and line.startswith("ATOM") and line[16] in " A":
            cur.append(line)
    return models


def coords(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])


def core_rmsf(models):
    """Per-residue CA RMSF after iterated superposition on the better-ordered half."""
    X = np.array([[[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                   for l in m if l[12:16] == " CA "] for m in models])
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
    return rmsf


_BASES = {}


def base(code):
    if code not in _BASES:
        _BASES[code] = Base(code)
    return _BASES[code]


class Base:
    def __init__(self, code):
        self.code = code
        self.models = read_models(code)
        self.lines = self.models[0]
        self.names = [l[12:16].strip() for l in self.lines]
        keys = [(l[21], l[22:27]) for l in self.lines]
        order = list(dict.fromkeys(keys))
        self.ri = np.array([order.index(k) for k in keys])
        self.n = len(order)
        self.resn = {self.ri[i]: self.lines[i][17:20] for i in range(len(self.lines))}
        self.atom = {(self.ri[i], self.names[i]): i for i in range(len(self.lines))}
        self.nums = [int(l[22:26]) for l in self.lines if l[12:16] == " CA "]
        self.rmsf = core_rmsf(self.models)

    def chi_moves(self, k):
        """(mask, axis atom b, axis atom c) for every side-chain dihedral (a, b, c, d) of
        residue k: the atoms on the c side of the b-c bond, found in the covalent graph of
        the residue (bonds < 1.95 A, < 1.3 A with a hydrogen); ring torsions are skipped."""
        from resuper.core import _CHI
        idx = [i for i in range(len(self.lines)) if self.ri[i] == k]
        X = coords([self.lines[i] for i in idx])
        names = [self.names[i] for i in idx]
        isH = [n.startswith("H") or (n[:1].isdigit() and n[1:2] == "H") for n in names]
        nb = {j: set() for j in range(len(idx))}
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                lim = 1.3 if isH[a] or isH[b] else 1.95
                if np.linalg.norm(X[a] - X[b]) < lim:
                    nb[a].add(b)
                    nb[b].add(a)
        pos = {n: j for j, n in enumerate(names)}
        out = []
        for quad in _CHI.get(self.resn[k], ()):
            if not all(q in pos for q in quad):
                continue
            b, c = pos[quad[1]], pos[quad[2]]
            seen, stack = {c}, [c]
            while stack:
                for y in nb[stack.pop()]:
                    if y != b and y not in seen:
                        seen.add(y)
                        stack.append(y)
            if any(names[j] in ("N", "CA", "C", "O") for j in seen):
                continue                    # ring through the backbone (proline)
            mask = np.zeros(len(self.lines), dtype=bool)
            mask[[idx[j] for j in seen if j != c]] = True
            out.append((mask, idx[b], idx[c]))
        return out

    def torsion_moves(self, k):
        """(mask, axis atom a, axis atom b) for phi_k (unless proline) and psi_k."""
        out = []
        if self.resn[k] != "PRO":
            m = (self.ri > k) | ((self.ri == k) &
                                 np.array([n not in H.N_TERM_H for n in self.names]))
            out.append((m, self.atom[(k, "N")], self.atom[(k, "CA")]))
        m = (self.ri > k) | ((self.ri == k) & np.array([n in H.O_ATOMS for n in self.names]))
        out.append((m, self.atom[(k, "CA")], self.atom[(k, "C")]))
        return out


# ---------------------------------------------------------------- design

def design():
    """{name: case}; case = dict(code, suite, scen, base, nconf, hinges [(k, sd)],
    states (k, [theta per conformer]) or None, flex [k], exp_flex bool)."""
    cases = {}
    ref_single = {}
    for code in proteins():
        p = f"{REF}/base_{code}.out"
        if os.path.exists(p):
            ref_single[code] = len(parse_report(open(p).read())) == 1
    for p, code in enumerate(proteins()):
        rng = np.random.default_rng(7000 + SEED_OFFSET + p)
        n = sum(1 for l in read_models(code)[0] if l[12:16] == " CA ")
        pos = lambda lo, hi: int(round(rng.uniform(lo, hi) * n))
        amp = lambda: float(rng.uniform(10, 20))
        frac = lambda lo, hi: int(round(rng.uniform(lo, hi) * n))
        c = lambda scen, suite, **kw: cases.__setitem__(
            f"{code}_{scen}", dict(code=code, suite=suite, scen=scen, base="syn", nconf=20,
                                   hinges=[], states=None, flex=[], exp_flex=False) | kw)
        c("c1", "clusters")
        c("c2", "clusters", hinges=[(pos(0.3, 0.7), amp())])
        c("c3", "clusters", hinges=[(pos(0.2, 0.4), amp()), (pos(0.6, 0.8), amp())])
        c("c4", "clusters", hinges=[(pos(0.15, 0.3), amp()), (pos(0.4, 0.6), amp()),
                                    (pos(0.7, 0.85), amp())])
        e2 = [(pos(0.3, 0.7), amp())]
        e3 = [(pos(0.2, 0.4), amp()), (pos(0.6, 0.8), amp())]
        if ref_single.get(code):
            for scen, hs in (("e1", []), ("e2", e2), ("e3", e3)):
                c(scen, "expclus", base="exp", hinges=hs, exp_flex=True)
        L = frac(0.25, 0.4)
        nterm = bool(rng.integers(0, 2))
        c("idr_term", "idr", flex=list(range(L)) if nterm else list(range(n - L, n)))
        a, b = frac(0.15, 0.25), frac(0.15, 0.25)
        c("idr_both", "idr", flex=list(range(a)) + list(range(n - b, n)))
        s, w = pos(0.35, 0.5), int(rng.integers(15, 26))
        c("idr_loop", "idr", flex=list(range(s, min(n, s + w))))
        c("idp_full", "idp", flex=list(range(n)))
        s1 = pos(0.3, 0.6)
        c("idp_seg", "idp", flex=[i for i in range(n) if not s1 <= i < s1 + 12])
        t1, t2 = pos(0.15, 0.3), pos(0.6, 0.75)
        c("idp_2seg", "idp", flex=[i for i in range(n)
                                   if not (t1 <= i < t1 + 12 or t2 <= i < t2 + 12)])
        hk = pos(0.35, 0.65)
        st = lambda counts, thetas: [t for cnt, t in zip(counts, thetas) for _ in range(cnt)]
        c("ms2_15", "multi", states=(hk, st([10, 10], [0, 15])))
        c("ms2_30", "multi", states=(hk, st([10, 10], [0, 30])))
        c("ms2_skew", "multi", states=(hk, st([16, 4], [0, 30])))
        c("ms3", "multi", states=(hk, st([7, 7, 6], [0, 25, -25])))
        h = (pos(0.35, 0.65), amp())
        for nc in (5, 10, 40):
            c(f"nc{nc}", "multi", hinges=[h], nconf=nc)
    return cases


def build_case(name, case, seed):
    B = base(case["code"])
    rng = np.random.default_rng(seed)
    base0 = coords(B.lines)
    sigma = np.clip(0.3 * B.rmsf, 0.05, 0.5)
    w = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    w /= np.sqrt((w ** 2).sum())
    moves = [(mv, sd) for k, sd in case["hinges"] for mv in B.torsion_moves(k)]
    moves += [(mv, IDR_SD) for k in case["flex"] for mv in B.torsion_moves(k)]
    moves += [(mv, CHI_SD) for k in case["flex"] for mv in B.chi_moves(k)]
    smoves = B.torsion_moves(case["states"][0]) if case["states"] else []
    nconf = case["nconf"] if case["base"] == "syn" else min(20, len(B.models))
    with open(f"{IN}/{name}.pdb", "w") as f:
        f.write(f"TITLE     benchmark B {name}\n")
        for cnum in range(nconf):
            X = base0.copy() if case["base"] == "syn" else coords(B.models[cnum])
            for (m, a, b), sd in moves:
                H.rotate(X, m, a, b, rng.normal(0, sd))
            if case["states"]:
                theta = case["states"][1][cnum]
                for m, a, b in smoves:
                    H.rotate(X, m, a, b, theta + rng.normal(0, 3.0))
            if case["base"] == "syn":
                X += rng.normal(0, 0.05, X.shape)
                z = rng.normal(0, 1, (B.n + 12, 3))
                field = np.stack([np.convolve(z[:, j], w, mode="valid") for j in range(3)], 1)
                X += (field * sigma[:, None])[B.ri]
            f.write(f"MODEL     {cnum + 1:4d}\n")
            for l, p in zip(B.lines, X):
                f.write(l[:16] + " " + l[17:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:])
            f.write("ENDMDL\n")
        f.write("END\n")


def truth(case):
    """{residue number: domain id | -1 flexible | None not scored}, number of domains."""
    B = base(case["code"])
    disordered = {i for i, v in enumerate(B.rmsf) if v > FLEX_RMSF} if case["exp_flex"] else set()
    hinge = {k for k, _ in case["hinges"]} | ({case["states"][0]} if case["states"] else set())
    flex = set(case["flex"])
    out, d, runs, cur = {}, 0, [], []
    for i in range(B.n):
        if i in hinge or i in flex:
            if cur:
                runs.append(cur)
            cur = []
        else:
            cur.append(i)
    if cur:
        runs.append(cur)
    for run in runs:
        rigid = [i for i in run if i not in disordered]
        lab = d if len(rigid) >= MU else -1
        d += lab >= 0
        for i in run:
            out[B.nums[i]] = lab if i not in disordered else -1
    for i in flex:
        out[B.nums[i]] = -1
    for i in hinge:
        out[B.nums[i]] = None
    return out, d


# ---------------------------------------------------------------- defective inputs

def split_file(path):
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


def write_models(path, head, models, tail, numbers=None, bare=False, crlf=False):
    out = list(head)
    for i, body in enumerate(models):
        if not bare:
            out.append(f"MODEL     {(numbers[i] if numbers else i + 1):4d}\n")
        out += body
        if not bare:
            out.append("ENDMDL\n")
    out += tail or ["END\n"]
    text = "".join(out)
    if crlf:
        text = text.replace("\n", "\r\n")
    open(path, "w", newline="").write(text)


def defect(code, kind):
    head, models, tail = split_file(f"{PDB}/{code}.pdb")
    models = [[l for l in m if l.startswith(("ATOM", "HETATM", "TER"))] for m in models]
    atoms = [l for l in models[0] if l.startswith("ATOM")]
    resnums = list(dict.fromkeys(int(l[22:26]) for l in atoms))
    mid = resnums[len(resnums) // 2]
    is_res = lambda l, r: l.startswith("ATOM") and int(l[22:26]) == r
    path = f"{IN}/err_{code}_{kind}.pdb"
    kw = {}
    if kind == "miss_sidechain_atom":
        r = next(r for r in resnums[len(resnums) // 2:]
                 if any(is_res(l, r) and l[12:16] == " CB " for l in models[0]))
        models[2] = [l for l in models[2] if not (is_res(l, r) and l[12:16] == " CB ")]
    elif kind == "miss_ca_all":
        models = [[l for l in m if not (is_res(l, mid) and l[12:16] == " CA ")] for m in models]
    elif kind == "miss_ca_one_model":
        models[2] = [l for l in models[2] if not (is_res(l, mid) and l[12:16] == " CA ")]
    elif kind == "gap3":
        models = [[l for l in m if not any(is_res(l, r) for r in (mid, mid + 1, mid + 2))]
                  for m in models]
    elif kind == "hetatm_ligand":
        def lig(m):
            x = [l for l in m if l.startswith("ATOM")][-1]
            extra = []
            for j in range(6):
                extra.append(f"HETATM{9000 + j:5d}  C{j}  LIG A 901    "
                             f"{float(x[30:38]) + j:8.3f}{float(x[38:46]):8.3f}"
                             f"{float(x[46:54]):8.3f}  1.00  0.00           C\n")
            extra.append(f"HETATM{9010:5d}  O   HOH A 902    {float(x[30:38]) + 9:8.3f}"
                         f"{float(x[38:46]):8.3f}{float(x[46:54]):8.3f}  1.00  0.00           O\n")
            return m + extra
        models = [lig(m) for m in models]
    elif kind == "mse":
        r = next((r for r in resnums if any(is_res(l, r) and l[17:20] == "MET" for l in atoms)),
                 resnums[1])
        def mse(l):
            if is_res(l, r):
                l = "HETATM" + l[6:17] + "MSE" + l[20:]
                if l[12:16] == " SD ":
                    l = l[:12] + "SE  " + l[16:76] + "SE" + l[78:]
            return l
        models = [[mse(l) for l in m] for m in models]
    elif kind == "altloc":
        def alt(m):
            out = []
            for l in m:
                if is_res(l, mid):
                    out.append(l[:16] + "A" + l[17:54] + "  0.60" + l[60:])
                    b = l[:16] + "B" + l[17:30] + f"{float(l[30:38]) + 0.3:8.3f}" + l[38:54] \
                        + "  0.40" + l[60:]
                    out.append(b)
                else:
                    out.append(l)
            return out
        models = [alt(m) for m in models]
    elif kind == "icode":
        def ic(l):
            if l.startswith("ATOM") and int(l[22:26]) in (mid + 1, mid + 2, mid + 3):
                code_ = "ABC"[int(l[22:26]) - mid - 1]
                return l[:22] + f"{mid:4d}" + code_ + l[27:]
            return l
        models = [[ic(l) for l in m] for m in models]
    elif kind == "strip_h":
        models = [[l for l in m if not (l.startswith("ATOM") and
                                        (l[76:78].strip() == "H" or l[12:16].strip()[:1] == "H"))]
                  for m in models]
    elif kind == "split_chain":
        def sc(m):
            out = [l[:21] + "B" + l[22:] if l.startswith("ATOM") and int(l[22:26]) > mid else l
                   for l in m if not l.startswith("TER")]
            k = next(i for i, l in enumerate(out) if l.startswith("ATOM") and l[21] == "B")
            return out[:k] + ["TER\n"] + out[k:] + ["TER\n"]
        models = [sc(m) for m in models]
    elif kind == "dup_atom":
        k = next(i for i, l in enumerate(models[1]) if is_res(l, mid) and l[12:16] == " CA ")
        models[1].insert(k, models[1][k])
    elif kind == "miss_residue_one_model":
        models[3] = [l for l in models[3] if not is_res(l, mid)]
    elif kind == "shuffled_atoms":
        def shuf(m):
            out, cur, key = [], [], None
            for l in m + [None]:
                k = (l[21], l[22:27]) if l and l.startswith("ATOM") else None
                if k != key and cur:
                    out += cur[::-1]
                    cur = []
                if l is None:
                    break
                if k:
                    cur.append(l)
                else:
                    out.append(l)
                key = k
            return out
        models[1] = shuf(models[1])
    elif kind == "crlf":
        kw["crlf"] = True
    elif kind == "model_numbers":
        kw["numbers"] = [5 + 3 * i for i in range(len(models))]
    elif kind == "far_model":
        models[1] = [l[:30] + f"{float(l[30:38]) + 900:8.3f}{float(l[38:46]) - 700:8.3f}"
                     f"{float(l[46:54]) + 500:8.3f}" + l[54:]
                     if l.startswith(("ATOM", "HETATM")) else l for l in models[1]]
    elif kind == "two_models":
        models = models[:2]
    elif kind == "single_model":
        models = models[:1]
    elif kind == "bare_single_model":
        models = models[:1]
        kw["bare"] = True
    elif kind == "identical_models":
        models = [models[0]] * len(models)
    elif kind == "nan_coord":
        k = next(i for i, l in enumerate(models[1]) if is_res(l, mid) and l[12:16] == " CA ")
        models[1][k] = models[1][k][:30] + "     nan" + models[1][k][38:]
    elif kind == "junk_lines":
        k = len(models[1]) // 2
        models[1] = models[1][:k] + ["THIS IS NOT A PDB RECORD\n", models[1][k][:40] + "\n"] \
            + models[1][k + 1:]
    elif kind == "empty":
        open(path, "w").close()
        return
    elif kind == "no_model_records":
        kw["bare"] = True
    elif kind == "tiny":
        keep = set(resnums[:6])
        models = [[l for l in m if l.startswith("ATOM") and int(l[22:26]) in keep]
                  for m in models]
    write_models(path, head, models, tail, **kw)


# ---------------------------------------------------------------- running

def all_inputs():
    cases = design()
    names = list(cases)
    for code in proteins():
        names.append(f"base_{code}")
    for code in proteins()[:8]:
        names += [f"err_{code}_{k}" for k in DEFECTS]
    return cases, names


RUNNER = ("import json, sys; import resuper.core as C; "
          "[setattr(C, k, v) for k, v in json.loads(sys.argv[1]).items()]; "
          "from resuper.cli import main; sys.exit(main(sys.argv[2:]))")


def run_one(job):
    """job = (method, name, config dict, output tag): resuper with module constants of
    resuper.core overridden, output in OUT/<method><tag>/."""
    method, name, cfg, tag = job
    d = f"{OUT}/{method}{tag}"
    os.makedirs(d, exist_ok=True)
    try:
        p = subprocess.run([sys.executable, "-c", RUNNER, json.dumps(cfg), "--method", method,
                            f"{IN}/{name}.pdb"], capture_output=True, text=True,
                           timeout=TIMEOUT, env=dict(os.environ, PYTHONPATH="src"))
        text, rc = p.stdout + p.stderr, p.returncode
    except subprocess.TimeoutExpired:
        text, rc = "TIMEOUT\n", 124
    open(f"{d}/{name}.out", "w").write(text)
    open(f"{d}/{name}.rc", "w").write(f"{rc}\n")


def load(method, name):
    if method.startswith("default") and TAG:
        method += TAG
    if method == "baseline":
        text = open(f"{REF}/{name}.out", errors="replace").read()
        rc = int(open(f"{REF}/{name}.rc").read().split()[0])
    else:
        text = open(f"{OUT}/{method}/{name}.out", errors="replace").read()
        rc = int(open(f"{OUT}/{method}/{name}.rc").read().split()[0])
    return text, rc


# ---------------------------------------------------------------- scoring

def ari(a, b):
    """Adjusted Rand index of two labelings (lists of hashables)."""
    from collections import Counter
    n = len(a)
    if n < 2:
        return 1.0
    comb = lambda x: x * (x - 1) / 2
    sab = sum(comb(v) for v in Counter(zip(a, b)).values())
    sa = sum(comb(v) for v in Counter(a).values())
    sb = sum(comb(v) for v in Counter(b).values())
    exp = sa * sb / comb(n)
    mx = (sa + sb) / 2
    return 1.0 if mx == exp else (sab - exp) / (mx - exp)


def score_case(labels, ndom, pred):
    s = H.score(labels, ndom, pred)
    scored = [r for r, l in labels.items() if l is not None]
    owner = {}
    for j, p in enumerate(pred):
        for r in p:
            owner.setdefault(r, j)
    s["ari"] = ari([labels[r] for r in scored], [owner.get(r, -1) for r in scored])
    s["selected"] = len(set().union(*pred) & set(scored)) / max(1, len(scored)) if pred else 0.0
    return s


def outcome(text, rc, clean_sets, expect):
    if rc in (124, 137, 139, -11, -9) or "Traceback" in text or "egmentation" in text \
            or "TIMEOUT" in text or rc < 0:
        return "crash"
    sets = [d[0] for d in parse_report(text)]
    if rc != 0 or (not sets and "ERROR" in text.upper()):
        return "refused"     # an error message (CYRANGE reports FATAL ERROR with exit 0)
    if not sets:
        return "empty"
    if expect == "any":
        return "answer"
    num = lambda ss: {i for _, i in set().union(*ss)} if ss else set()
    a, b = num(sets), num(clean_sets)
    j = len(a & b) / max(1, len(a | b))
    return "near" if len(sets) == len(clean_sets) and j >= 0.8 else "drift"


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "score"
    os.makedirs(IN, exist_ok=True)
    if cmd == "build-base":
        for code in proteins():
            head, models, tail = split_file(f"{PDB}/{code}.pdb")
            write_models(f"{IN}/base_{code}.pdb", head, models, tail)
    elif cmd == "build":
        cases, names = all_inputs()
        for s, (name, case) in enumerate(cases.items()):
            build_case(name, case, 31000 + 10 * SEED_OFFSET + s)
        for code in proteins():
            head, models, tail = split_file(f"{PDB}/{code}.pdb")
            write_models(f"{IN}/base_{code}.pdb", head, models, tail)
        for code in proteins()[:8]:
            for k in DEFECTS:
                defect(code, k)
        print(f"{len(names)} inputs in {IN}/")
    elif cmd == "capture":
        names = [a for a in sys.argv[2:] if a not in ("--bench", BSET)] or all_inputs()[1]
        os.makedirs(REF, exist_ok=True)
        subprocess.run(["ssh", H.ORACLE, f"rm -rf {REMOTE} && mkdir -p {REMOTE}"], check=True)
        subprocess.run(["scp", "-q"] + [f"{IN}/{n}.pdb" for n in names]
                       + [f"{H.ORACLE}:{REMOTE}/"], check=True)
        subprocess.run(["ssh", H.ORACLE, f"cd {REMOTE} && ls *.pdb | xargs -P 32 -I{{}} sh -c "
                        f"'b=$(basename {{}} .pdb); timeout {TIMEOUT} /usr/local/bin/cyrange {{}} "
                        "> $b.out 2>&1; echo $? > $b.rc'"], check=True)
        tar = subprocess.run(["ssh", H.ORACLE, f"cd {REMOTE} && tar cf - *.out *.rc"],
                             capture_output=True, check=True).stdout
        subprocess.run(["tar", "xf", "-", "-C", REF], input=tar, check=True)
        print(f"{len(names)} reference reports in {REF}/")
    elif cmd == "run":
        names = all_inputs()[1]
        args = [a for a in sys.argv[2:] if "=" not in a and a != "--bench" and a != BSET]
        cfg = {a.split("=", 1)[0]: ast.literal_eval(a.split("=", 1)[1])
               for a in sys.argv[2:] if "=" in a and not a.startswith("tag=")}
        tag = next((a[4:] for a in sys.argv[2:] if a.startswith("tag=")), "")
        if args:                # only inputs whose name starts with one of the prefixes
            names = [n for n in names if n.startswith(tuple(args))]
        methods = ("default",) if cfg or tag else ("reference", "default")
        jobs = [(m, n, cfg, tag) for m in methods for n in names]
        with Pool() as pool:
            pool.map(run_one, jobs, chunksize=4)
        print(f"{len(jobs)} resuper runs in {OUT}/")
    else:
        global TAG
        TAG = next((a[4:] for a in sys.argv[2:] if a.startswith("tag=")), "")
        only = next((a[7:].split(",") for a in sys.argv[2:] if a.startswith("suites=")), None)
        score(*all_inputs(), only)


def score(cases, names, only=None):
    """only: suite names to score ("errors" for the defective inputs); None = all."""
    if only:
        cases = {n: c for n, c in cases.items() if c["suite"] in only}
    print(f"proteins: {len(proteins())}; ground-truth ensembles: {len(cases)}")
    rows = {}
    for name, case in cases.items():
        labels, ndom = truth(case)
        rows[name] = (case, ndom, {m: score_case(labels, ndom, [{i for _, i in d[0]} for d in
                                   parse_report(load(m, name)[0])]) for m in METHODS})
    hdr = "count  mix   ARI  contam  Jacc"
    print(f"\n{'scenario':10s} {'n':>3s}  " + "  ".join(f"{m:>30s}" for m in METHODS))
    print(" " * 16 + "  ".join(f"{hdr:>30s}" for _ in METHODS))
    for suite, scens in SUITES.items():
        groups = [(s, [n for n, (c, _, _) in rows.items() if c["scen"] == s]) for s in scens]
        groups.append((f"[{suite}]", [n for n, (c, _, _) in rows.items() if c["suite"] == suite]))
        for g, ns in groups:
            if not ns:
                continue
            cells = []
            for m in METHODS:
                r = [rows[n][2][m] for n in ns]
                cells.append(f"{sum(x['count'] for x in r):4d} {sum(x['mixing'] for x in r):4d} "
                             f"{np.mean([x['ari'] for x in r]):5.2f}  "
                             f"{np.mean([x['contam'] for x in r]):5.3f} "
                             f"{np.mean([x['jaccard'] for x in r]):5.3f}")
            print(f"{g:10s} {len(ns):3d}  " + "  ".join(f"{c:>30s}" for c in cells))
    if only and "errors" not in only:
        if "idp" not in only:
            return
    ns = [n for n, (c, _, _) in rows.items() if c["scen"] == "idp_full"]
    print("\nidp_full (no true domain): mean domains reported / mean fraction of residues "
          "selected / ensembles with no domain")
    for m in METHODS:
        r = [rows[n][2][m] for n in ns]
        print(f"  {m:10s} {np.mean([x['n'] for x in r]):.2f} / "
              f"{np.mean([x['selected'] for x in r]):.2f} / {sum(x['n'] == 0 for x in r)}")
    if only and "errors" not in only:
        return
    print("\nerrors: outcome counts per defect (8 proteins); expectation in brackets")
    print(f"{'defect':24s}" + "".join(f"{m:>34s}" for m in METHODS))
    tot = {m: {} for m in METHODS}
    for k, expect in DEFECTS.items():
        cells = []
        for m in METHODS:
            cnt = {}
            for code in proteins()[:8]:
                clean = [d[0] for d in parse_report(load(m, f"err_{code}_clean")[0])]
                o = outcome(*load(m, f"err_{code}_{k}"), clean, expect)
                cnt[o] = cnt.get(o, 0) + 1
                tot[m][o] = tot[m].get(o, 0) + 1
            cells.append(" ".join(f"{o}:{c}" for o, c in sorted(cnt.items())))
        print(f"{k:18s}[{expect:4s}]" + "".join(f"{c:>34s}" for c in cells))
    print("totals: " + " | ".join(f"{m}: " + " ".join(f"{o}:{c}" for o, c in sorted(t.items()))
                                  for m, t in tot.items()))


if __name__ == "__main__":
    main()
