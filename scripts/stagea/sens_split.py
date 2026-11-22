"""Sensitivity of the sound Stage A to the split-test thresholds and the tree linkage.

For every evaluation case (real proteins with reference reports, plus the synthetic
intact-hinge ensembles with oracle default reports in cywork/seg/) this builds the core
set (published Q_i cutoff), the clustering tree for each linkage, and caches every
cluster RMSD; the split walk of cyrange._rmsd_split_groups is then replayed for a grid
of (SPLIT_TAU, SPLIT_REL). Reported per setting: how many cases get the reference's
domain count (Stage-A groups; Stage B's drop of domains under mu residues is ignored).

usage: /usr/bin/python3 scripts/stagea/sens_split.py
"""
import glob, os, re, sys
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.abspath("scripts"))
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from eval_cases import CURATED, AISAR, pacsy_ids
from resuper.core import (order_parameters, _published_s_cut, _ca_variance_matrix,
                             Params)
from resuper.geometry import rmsd_to_mean
from resuper.pdb import parse_pdb

TAUS = [0.06, 0.08, 0.10, 0.11, 0.12, 0.14, 0.16, 0.20]
RELS = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7]
LINKAGES = ["ward", "average_sd"]
MU = Params().min_cluster


def cases():
    for base in dict.fromkeys(CURATED + AISAR + pacsy_ids()):
        out, pdb = f"testdata/{base}.out", f"testdata/{base}.pdb"
        if os.path.exists(out) and os.path.exists(pdb) and "Optimal range" in open(out).read():
            yield "real", base, pdb, len(re.findall("Optimal range", open(out).read()))
    for pat in ("ht_*", "hz_*", "hu_*", "hn_*", "hp_*", "hg_*", "fh_*"):
        for f in sorted(glob.glob(f"cywork/seg/{pat}.default")):
            yield "syn", os.path.basename(f)[:-8], f[:-8] + ".pdb", \
                len(re.findall("Optimal range", open(f).read()))


def tree(V, kind):
    if kind == "average_sd":
        D = np.sqrt(np.maximum(V, 0.0))
        np.fill_diagonal(D, 0.0)
        return to_tree(linkage(squareform((D + D.T) / 2, checks=False), method="average"))
    return to_tree(linkage(squareform(V, checks=False), method="ward"))


def groups(root, rm, tau, rel):
    out = []

    def walk(node, parent_m):
        ids = tuple(sorted(node.pre_order()))
        if node.is_leaf() or len(ids) < MU:
            if len(ids) >= MU:
                out.append(ids)
            return
        L = tuple(sorted(node.left.pre_order()))
        R = tuple(sorted(node.right.pre_order()))
        big = [c for c, g in ((node.left, L), (node.right, R)) if len(g) >= MU]
        if len(big) == 1:
            walk(big[0], parent_m)
            return
        m = rm(ids) - rm(L) - rm(R)
        if len(big) == 2:
            if m >= tau and (parent_m is None or m >= rel * parent_m):
                walk(node.left, m)
                walk(node.right, m)
            else:
                out.append(ids)
        elif m < tau:
            out.append(ids)
    walk(root, None)
    return out


def main():
    score = {}
    total = {"real": 0, "syn": 0}
    for kind, name, pdb, want in cases():
        st = parse_pdb(pdb)
        order = order_parameters(st)
        cut = _published_s_cut(order)
        core = [i for i, d in enumerate(order) if any(v >= cut for v in d.values())]
        total[kind] += 1
        if len(core) < MU:
            for key in [(l, t, r) for l in LINKAGES for t in TAUS for r in RELS]:
                score.setdefault((kind,) + key, 0)
                score[(kind,) + key] += want == 1
            continue
        V = _ca_variance_matrix(st, core)
        res = st.residues
        X = np.array([[[c[(res[i].chain, res[i].num)][a] for a in ("N", "CA", "C")]
                       for i in core] for c in st.conformers])
        cache = {}

        def rm(g):
            if g not in cache:
                e = X[:, list(g)].reshape(X.shape[0], -1, 3)
                cache[g] = rmsd_to_mean(e, ref_seq=e[0], iters=1)
            return cache[g]
        for lk in LINKAGES:
            root = tree(V, lk) if V.sum() > 1e-12 else None
            for t in TAUS:
                for r in RELS:
                    n = len(groups(root, rm, t, r)) if root is not None else 1
                    n = max(n, 1)
                    score.setdefault((kind, lk, t, r), 0)
                    score[(kind, lk, t, r)] += n == want
    for lk in LINKAGES:
        print(f"\n== linkage {lk}: domain count = reference (real /{total['real']}, "
              f"synthetic /{total['syn']})")
        print("tau \\ rel " + " ".join(f"{r:>11}" for r in RELS))
        for t in TAUS:
            print(f"{t:<9} " + " ".join(f"{score[('real', lk, t, r)]:>4}/{score[('syn', lk, t, r)]:<3}   "
                                         for r in RELS))


if __name__ == "__main__":
    main()
