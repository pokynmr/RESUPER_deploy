"""Relative motion of RESUPER segments that lie inside one deposited HELIX record (S16).

For every entry of results/expanded/evidence.csv with helix_split = True, and every pair of
adjacent RESUPER ranges inside one HELIX record:
  scatter   superpose each conformer on the first segment (Calpha, to conformer 1) and take the
            mean angle between the second segment's axis and its mean direction (degrees);
  control   the same quantity for the two halves of the longer segment (axis noise of short
            pieces of the same ensemble);
  hbond     lowest fraction of conformers keeping an O(i)-N(i+4) <= 3.5 A hydrogen bond across
            the junction (helix continuity);
  bend      angle between the two segment axes, mean and SD over conformers.

usage (from benchmarks/cyrange_comparison): python3 selection_study/scripts/helix_split_motion.py
PDB files are downloaded from the RCSB PDB into data/raw_pdb/ when missing.
Output: selection_study/results/expanded/helix_split_motion.csv
"""
import csv, os, re, sys, urllib.request
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402
sys.path.insert(0, cfg.REPO_SRC)
from resuper.pdb import parse_pdb  # noqa: E402
from resuper.geometry import _kabsch  # noqa: E402

RES = os.path.join(cfg.BENCH_DIR, "selection_study", "results", "expanded")
OUT = os.path.join(RES, "helix_split_motion.csv")


def helices(path):
    return [(int(l[21:25]), int(l[33:37])) for l in open(path) if l.startswith("HELIX ")]


def segments(s):
    return sorted((int(a), int(b)) for part in s.split("|")
                  for a, b in re.findall(r"(-?\d+)\.\.(-?\d+)", part.split("(")[0]))


def axis(P):
    u = np.linalg.svd(P - P.mean(0))[2][0]
    return u if np.dot(u, P[-1] - P[0]) > 0 else -u


def scatter(A, B):
    ax = []
    for k in range(len(A)):
        U = _kabsch(A[k] - A[k].mean(0), A[0] - A[0].mean(0))
        ax.append(axis((B[k] - A[k].mean(0)) @ U.T))
    ax = np.array(ax); m = ax.mean(0); m /= np.linalg.norm(m)
    return float(np.degrees(np.arccos(np.clip(ax @ m, -1, 1))).mean())


def main():
    os.makedirs(cfg.RAW_PDB_DIR, exist_ok=True)
    rows = [r for r in csv.DictReader(open(os.path.join(RES, "evidence.csv"))) if r["helix_split"] == "True"]
    out = []
    for r in rows:
        pid = r["pdb_id"]; path = os.path.join(cfg.RAW_PDB_DIR, f"{pid}.pdb")
        if not os.path.exists(path):
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.pdb", path)
        st = parse_pdb(path); C = st.conformers; ch = st.residues[0].chain
        ok = lambda n: (ch, n) in C[0]
        ca = lambda ns: np.array([[c[(ch, n)]["CA"] for n in ns] for c in C], float)
        H = helices(path); seg = segments(r["ranges_resuper"])
        for (a1, b1), (a2, b2) in zip(seg, seg[1:]):
            hh = [h for h in H if h[0] - 3 <= a1 and a2 <= h[1] - 4]
            if not hh:
                continue
            L = [n for n in range(a1, b1 + 1) if ok(n)]; R = [n for n in range(a2, b2 + 1) if ok(n)]
            A, B = ca(L), ca(R)
            bend = np.array([np.degrees(np.arccos(np.clip(np.dot(axis(x), axis(y)), -1, 1))) for x, y in zip(A, B)])
            lg = A if len(L) >= len(R) else B; h = lg.shape[1] // 2
            ctl = scatter(lg[:, :h], lg[:, h:])
            hb = []
            for i in range(b1 - 4, a2 + 1):
                if ok(i) and ok(i + 4) and "O" in C[0][(ch, i)] and "N" in C[0][(ch, i + 4)]:
                    hb.append(np.mean([np.linalg.norm(np.subtract(c[(ch, i)]["O"], c[(ch, i + 4)]["N"])) <= 3.5 for c in C]))
            sc = scatter(A, B)
            out.append(dict(pdb_id=pid, helix=f"{hh[0][0]}-{hh[0][1]}", segment_1=f"{a1}-{b1}", segment_2=f"{a2}-{b2}",
                            n_models=len(C), scatter_deg=round(sc, 1), control_deg=round(ctl, 1),
                            ratio=round(sc / ctl, 2), junction_hbond_min=round(min(hb), 2) if hb else "",
                            bend_mean_deg=round(float(bend.mean()), 1), bend_sd_deg=round(float(bend.std()), 1)))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print(f"{len(rows)} entries, {len(out)} junctions; scatter > control in "
          f"{sum(o['scatter_deg'] > o['control_deg'] for o in out)}; wrote {OUT}")


if __name__ == "__main__":
    main()
