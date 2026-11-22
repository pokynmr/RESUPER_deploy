"""Phase 7: structural analysis of the proteins where the methods (or the paper) differ.

Questions
  1VVD / 1VVE  one rigid body or two that move relative to each other?
  2D21         does reference-mode domain 1 span both lobes?
  2KR6         why do both methods differ from the paper's CYRANGE values?
  1AEY, 1AIW, 1PKT  are the residues the methods disagree on flexible?

Analyses (RESUPER is only called or imported, never modified)
  RMSF       per-residue Cα fluctuation after superposing all models on a chosen
             set of residues (iterative fit to the mean). High RMSF = mobile.
  rotation   "one body or two": superpose every model on domain 1, then find the
             rotation that best superposes its domain 2 onto the mean domain 2.
             The angle of that rotation, per model, measures how much domain 2
             swings relative to domain 1. Backbone atoms N, CA, C.
  baseline   the same rotation angle between the two sequence halves of every
             single-domain protein on which both methods agree. Halves of one
             rigid domain show the angle produced by noise alone, so the
             two-domain angles are judged against real data, not a guessed cutoff.
  states     `resuper --states <ID>.pdb` (RESUPER's multiconformer detector).
  SS         helix / strand ranges from the HELIX / SHEET records of the PDB file.

Outputs: results/phase7_structure/
  interdomain.csv      one row per two-body pair (rotation angles, RMSDs)
  baseline_halves.csv  the baseline angles
  rmsf/<ID>.csv        per-residue RMSF (several frames)
  states.csv           state-detector results
  figures: fig5_rotation.png, fig6_rmsf_<ID>.png, fig7_ensemble_<ID>.png
  summary.txt

usage (from benchmarks/cyrange_comparison):  python scripts/structural_analysis.py
"""
import csv
import os
import subprocess
import sys
import tempfile
import shutil
import statistics as stt

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import Patch, Rectangle

import bench_config as cfg
from report_parser import parse_report

RUNS = os.path.join(cfg.RESULTS_DIR, "phase4_runs")
CMP = os.path.join(cfg.RESULTS_DIR, "phase5_comparison")
OUT = os.path.join(cfg.RESULTS_DIR, "phase7_structure")
A, B = "resuper_default", "reference_mode"
DOMAIN = ["#2a78d6", "#eb6834"]
OTHER, INK, INK2, MUTED = "#c9c8c2", "#0b0b0b", "#52514e", "#8a897f"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})

cfg.import_resuper()
from resuper.pdb import parse_pdb, BACKBONE          # noqa: E402
from resuper.geometry import superpose               # noqa: E402


# ---------------------------------------------------------------- geometry
def kabsch(P, Q):
    """Rotation R and translation t with P @ R.T + t ≈ Q (least squares)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - pc @ R.T


def angle_deg(R):
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1.0, 1.0))))


class Ensemble:
    def __init__(self, pdb_id):
        self.id = pdb_id
        self.st = parse_pdb(os.path.join(cfg.RAW_PDB_DIR, f"{pdb_id}.pdb"))
        self.keys = [(r.chain, r.num) for r in self.st.residues]
        self.index = {k: i for i, k in enumerate(self.keys)}
        K, n = self.st.n_models, len(self.keys)
        X = np.zeros((K, n, 3, 3))                   # model, residue, N/CA/C, xyz
        for k, conf in enumerate(self.st.conformers):
            for i, key in enumerate(self.keys):
                for j, a in enumerate(BACKBONE):
                    X[k, i, j] = conf[key][a]
        self.X = X
        chains = sorted({c for c, _ in self.keys})
        self.default_chain = chains[0] if len(chains) == 1 else ""

    def idx(self, residues):
        return sorted(self.index[r] for r in residues if r in self.index)

    def fit(self, fit_idx, iters=50):
        """All models superposed on backbone atoms of fit_idx, iterated to the mean."""
        K = self.X.shape[0]
        flat = self.X.reshape(K, -1, 3)
        sel = self.X[:, fit_idx].reshape(K, -1, 3)
        ref = sel[0]
        for _ in range(iters):
            moved = np.empty_like(flat)
            for k in range(K):
                R, t = kabsch(sel[k], ref)
                moved[k] = flat[k] @ R.T + t
            new_ref = moved.reshape(self.X.shape)[:, fit_idx].reshape(K, -1, 3).mean(axis=0)
            if np.abs(new_ref - ref).max() < 1e-8:
                break
            ref = new_ref
        return moved.reshape(self.X.shape)

    def rmsf_ca(self, fit_idx):
        Y = self.fit(fit_idx)[:, :, 1]                # CA
        return np.sqrt(((Y - Y.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))

    def rotation(self, idx1, idx2):
        """Per-model rotation angle (deg) and centroid shift (Å) of body 2 after
        superposing body 1."""
        Y = self.fit(idx1)
        D2 = Y[:, idx2].reshape(Y.shape[0], -1, 3)
        mean2 = D2.mean(axis=0)
        ang, shift = [], []
        for k in range(Y.shape[0]):
            R, _ = kabsch(D2[k], mean2)
            ang.append(angle_deg(R))
            shift.append(float(np.linalg.norm(D2[k].mean(axis=0) - mean2.mean(axis=0))))
        return np.array(ang), np.array(shift)


def check_kabsch():
    """Our superposition must agree with RESUPER's own superpose()."""
    rng = np.random.default_rng(0)
    P = rng.normal(size=(20, 3))
    th = 0.7
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    Q = P @ Rz.T + np.array([1.0, -2.0, 0.5]) + rng.normal(scale=0.05, size=P.shape)
    R, t = kabsch(P, Q)
    return float(np.abs((P @ R.T + t) - superpose(P, Q)).max())


# ---------------------------------------------------------------- inputs
def domains(pdb_id, method, ens):
    rep = parse_report(open(os.path.join(RUNS, pdb_id, f"{method}.out")).read())
    return rep, [d.residue_set(ens.default_chain) for d in rep.domains]


def ss_ranges(pdb_id):
    """Helix and strand residue ranges from HELIX / SHEET records."""
    out = []
    for line in open(os.path.join(cfg.RAW_PDB_DIR, f"{pdb_id}.pdb")):
        try:
            if line.startswith("HELIX "):
                out.append(("helix", line[19].strip() or "A", int(line[21:25]), int(line[33:37])))
            elif line.startswith("SHEET "):
                out.append(("strand", line[21].strip() or "A", int(line[22:26]), int(line[33:37])))
        except ValueError:
            pass
    return out


def run_states(pdb_id):
    with tempfile.TemporaryDirectory() as td:
        shutil.copyfile(os.path.join(cfg.RAW_PDB_DIR, f"{pdb_id}.pdb"), os.path.join(td, f"{pdb_id}.pdb"))
        env = dict(os.environ, PYTHONPATH=cfg.REPO_SRC, PYTHONDONTWRITEBYTECODE="1")
        p = subprocess.run([sys.executable, "-m", "resuper", "--states", f"{pdb_id}.pdb"],
                           cwd=td, env=env, capture_output=True, text=True, timeout=3600)
    line = next((l.strip() for l in p.stdout.splitlines() if "Conformational states" in l), "")
    return p.returncode, line


# ---------------------------------------------------------------- figures
def fig_rmsf(ens, frames, tracks, ss, diff_res, title, path):
    """RMSF profiles (one line per fitting frame) above domain tracks and SS."""
    n = len(ens.keys)
    x = np.arange(n)
    fig, (ax, axt) = plt.subplots(2, 1, figsize=(7.0, 3.9), dpi=200, sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.5], "hspace": 0.35})
    styles = ["-", "--", ":"]
    for (label, r), ls in zip(frames, styles):
        ax.plot(x, r, ls, color=INK if ls == "-" else INK2, lw=1.2, label=label)
    for i in diff_res:
        ax.axvspan(i - 0.5, i + 0.5, color="#fde2d6", lw=0, zorder=0)
    ax.set_yscale("log")
    ax.set_yticks([0.03, 0.1, 0.3, 1, 3, 10, 30, 100])
    ax.set_yticklabels(["0.03", "0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    lo = min(0.1, 0.8 * min(float(r.min()) for _, r in frames))   # never cut data off
    ax.set_ylim(lo, max(3.0, max(float(r.max()) for _, r in frames) * 1.3))
    ax.set_ylabel("Cα RMSF (Å, log scale)")
    ax.grid(axis="y", color="#ecebe6", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.03),
              ncol=3, handlelength=2.5)
    ax.set_title(title, loc="left", fontsize=8.5, color=INK, fontweight="bold")
    # tracks: SS, method A, method B
    labels = []
    for row, (name, per) in enumerate(tracks):
        y = len(tracks) + 1 - row          # method tracks above; SS row at y = 1
        for i, v in enumerate(per):
            axt.add_patch(Rectangle((i - 0.5, y - 0.35), 1, 0.7, lw=0,
                                    color=OTHER if v is None else (DOMAIN[v % 2] if v != "both" else "#1baf7a")))
        labels.append((y, name))
    for kind, ch, a, b in ss:
        ia = [i for i, k in enumerate(ens.keys) if k[0] == ch and a <= k[1] <= b]
        if ia:
            axt.add_patch(Rectangle((ia[0] - 0.5, 0.65), len(ia), 0.7, lw=0,
                                    color=INK if kind == "helix" else MUTED,
                                    hatch=None if kind == "helix" else "////", fill=kind == "helix"))
    labels.append((1, "helix (black) / strand (grey)"))
    axt.set_yticks([y for y, _ in labels])
    axt.set_yticklabels([l for _, l in labels], fontsize=6.5)
    axt.set_ylim(0.4, len(tracks) + 1.6)
    ticks = list(range(0, n, max(1, n // 10)))
    axt.set_xticks(ticks)
    axt.set_xticklabels([f"{ens.keys[t][1]}" for t in ticks], fontsize=6.5)
    axt.set_xlabel("Residue number" + ("  (shaded above: residues the methods assign differently)" if diff_res else ""))
    for s in ("right", "top", "left"):
        axt.spines[s].set_visible(False)
    axt.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.93, bottom=0.12)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def fig_ensemble(ens, bodies, names, title, path):
    """Cα traces of all models, superposed on each body in turn (one panel each)."""
    fig, axes = plt.subplots(1, len(bodies), figsize=(3.4 * len(bodies), 3.5), dpi=200)
    owner = {}
    for b, res in enumerate(bodies):
        for i in res:
            owner[i] = b
    for ax, (b, res) in zip(np.atleast_1d(axes), enumerate(bodies)):
        Y = ens.fit(res)[:, :, 1]
        mean = Y.mean(axis=0)
        c = mean.mean(axis=0)
        _, _, Vt = np.linalg.svd(mean - c)            # view: first two principal axes
        P = (Y - c) @ Vt[:2].T
        for k in range(P.shape[0]):
            for i in range(P.shape[1] - 1):
                col = DOMAIN[owner[i]] if owner.get(i) == owner.get(i + 1) and i in owner else OTHER
                ax.plot(P[k, i:i + 2, 0], P[k, i:i + 2, 1], color=col, lw=0.35, alpha=0.55)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"superposed on {names[b]}", fontsize=7.5, color=INK2)
    fig.suptitle(title, x=0.01, ha="left", fontsize=8.5, color=INK, fontweight="bold")
    fig.legend(handles=[Patch(color=DOMAIN[i], label=names[i]) for i in range(len(bodies))]
               + [Patch(color=OTHER, label="other residues")], loc="lower center",
               ncol=len(bodies) + 1, frameon=False, fontsize=6.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def fig_rotation(rows, base, path):
    fig, ax = plt.subplots(figsize=(7.0, 3.1), dpi=200)
    labels = [f"{r['pdb_id']} {r['bodies']}" for r in rows]
    for y, r in enumerate(rows):
        a = r["_angles"]
        ax.scatter(a, np.full(len(a), y) + np.random.default_rng(y).uniform(-0.15, 0.15, len(a)),
                   s=8, color=DOMAIN[0], alpha=0.6, lw=0)
        ax.plot([np.median(a)] * 2, [y - 0.3, y + 0.3], color=INK, lw=1.5)
    hi = np.percentile(base, 95)
    ax.axvspan(0, hi, color="#ecebe6", zorder=0)
    ax.text(hi + 2, -0.75, f"grey band: 95% of baseline angles\n(halves of single domains) are below {hi:.1f}°",
            fontsize=6.5, color=INK2, va="center")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_ylim(len(rows) - 0.4, -1.2)
    ax.set_xlabel("Rotation of body 2 relative to body 1 per model (degrees); bar = median")
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_title("Relative motion of the two bodies", loc="left",
                 fontsize=8.5, color=INK, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


# ---------------------------------------------------------------- main
def assignment(ens, sets):
    per = [None] * len(ens.keys)
    for j, s in enumerate(sets):
        for r in s:
            if r in ens.index:
                i = ens.index[r]
                per[i] = "both" if per[i] is not None else j
    return per


def main():
    log = cfg.setup_logging("structure")
    os.makedirs(os.path.join(OUT, "rmsf"), exist_ok=True)
    err = check_kabsch()
    log.info("superposition check against RESUPER's superpose(): max difference %.2e Å", err)
    assert err < 1e-6, "our superposition disagrees with RESUPER's"

    comp = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(CMP, "per_protein_comparison.csv")))}
    L = [f"superposition check vs RESUPER superpose(): max difference {err:.1e} Å"]

    # ---- baseline: halves of single domains where both methods agree on 1 domain
    base_rows = []
    for pid, c in comp.items():
        if not (c["ndom_a"] == "1" and c["ndom_b"] == "1"):
            continue
        ens = Ensemble(pid)
        _, sa = domains(pid, A, ens)
        idx = ens.idx(sa[0])
        if len(idx) < 16:
            continue
        h1, h2 = idx[: len(idx) // 2], idx[len(idx) // 2:]
        ang, _ = ens.rotation(h1, h2)
        base_rows.append({"pdb_id": pid, "n_residues": len(idx), "median_deg": round(float(np.median(ang)), 2),
                          "max_deg": round(float(ang.max()), 2), "_angles": ang})
    base_all = np.concatenate([r["_angles"] for r in base_rows])
    p95 = float(np.percentile(base_all, 95))
    L.append(f"\n== baseline: halves of {len(base_rows)} single-domain proteins "
             f"({base_all.size} model angles): median {np.median(base_all):.1f} deg, "
             f"95th percentile {p95:.1f} deg, max {base_all.max():.1f} deg")
    L.append("   per-protein median angles: " + ", ".join(
        f"{r['pdb_id']} {r['median_deg']}" for r in sorted(base_rows, key=lambda r: -r["median_deg"])))
    import json
    with open(os.path.join(OUT, "baseline_stats.json"), "w") as fh:
        json.dump({"n_proteins": len(base_rows), "n_angles": int(base_all.size),
                   "median_deg": round(float(np.median(base_all)), 1),
                   "p95_deg": round(p95, 1), "max_deg": round(float(base_all.max()), 1)}, fh, indent=1)
    with open(os.path.join(OUT, "baseline_halves.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb_id", "n_residues", "median_deg", "max_deg"],
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(base_rows)

    # ---- two-body pairs: every protein with 2 domains in either method
    pairs = []
    for pid, c in comp.items():
        for meth in (A, B):
            if c[f"ndom_{'a' if meth == A else 'b'}"] == "2":
                pairs.append((pid, meth))
    # 2D21 reference-mode domains overlap; the lobe question is answered with
    # RESUPER's (non-overlapping) domains, so skip the overlapping definition.
    pairs = [p for p in pairs if p != ("2D21", B)]
    seen, inter = set(), []
    for pid, meth in pairs:
        ens = Ensemble(pid)
        rep, sets = domains(pid, meth, ens)
        key = (pid, tuple(sorted(tuple(sorted(s)) for s in sets)))
        if key in seen:                 # both methods gave the same two domains
            continue
        seen.add(key)
        other = B if meth == A else A
        _, osets = domains(pid, other, ens)
        same = sorted(map(sorted, osets)) == sorted(map(sorted, sets))
        source = "both methods" if same else meth.replace("_", " ")
        i1, i2 = ens.idx(sets[0]), ens.idx(sets[1])
        ang, shift = ens.rotation(i1, i2)
        frac = float(np.mean(ang > p95))
        inter.append({
            "pdb_id": pid, "domains_from": source,
            "body1": rep.domains[0].text, "body2": rep.domains[1].text,
            "bodies": f"[{rep.domains[0].text}] vs [{rep.domains[1].text}]",
            "n_models": len(ang), "median_deg": round(float(np.median(ang)), 1),
            "p90_deg": round(float(np.percentile(ang, 90)), 1), "max_deg": round(float(ang.max()), 1),
            "median_shift_A": round(float(np.median(shift)), 2),
            "fraction_models_above_baseline_p95": round(frac, 2),
            "rmsd_body1": rep.domains[0].rmsd, "rmsd_body2": rep.domains[1].rmsd, "_angles": ang})
    inter.sort(key=lambda r: -r["median_deg"])
    with open(os.path.join(OUT, "interdomain.csv"), "w", newline="") as fh:
        f = [k for k in inter[0] if not k.startswith("_")]
        w = csv.DictWriter(fh, fieldnames=f, extrasaction="ignore")
        w.writeheader()
        w.writerows(inter)
    L.append("\n== two-body rotation (body 2 relative to body 1)")
    for r in inter:
        L.append(f"   {r['pdb_id']} ({r['domains_from']}): {r['bodies']}  median {r['median_deg']} deg, "
                 f"90th pct {r['p90_deg']}, max {r['max_deg']}; centroid shift median {r['median_shift_A']} A; "
                 f"{int(100 * r['fraction_models_above_baseline_p95'])}% of models above baseline 95th pct")
    fig_rotation(inter, base_all, os.path.join(OUT, "fig5_rotation.png"))

    # ---- case studies: RMSF + ensemble pictures
    cases = {
        "1VVD": dict(frames=["union"], bodies_from=B),
        "1VVE": dict(frames=["union"], bodies_from=A),
        "2D21": dict(frames=["union"], bodies_from=A),
        "2KR6": dict(frames=["union"], bodies_from=A),
        "1AEY": dict(frames=[], bodies_from=A),
        "1AIW": dict(frames=[], bodies_from=A),
        "1PKT": dict(frames=[], bodies_from=A),
    }
    L.append("\n== per-residue flexibility at the residues the methods assign differently")
    for pid, spec in cases.items():
        ens = Ensemble(pid)
        repA, sA = domains(pid, A, ens)
        repB, sB = domains(pid, B, ens)
        bodies_rep, bodies = (repA, sA) if spec["bodies_from"] == A else (repB, sB)
        per_a, per_b = assignment(ens, sA), assignment(ens, sB)
        diff = [i for i in range(len(ens.keys)) if per_a[i] != per_b[i]]
        frames = []
        for j, s in enumerate(bodies):
            frames.append((f"fit on {spec['bodies_from'].replace('_', ' ')} domain {j + 1}",
                           ens.rmsf_ca(ens.idx(s))))
        if len(bodies) > 1:
            u = set().union(*bodies)
            frames.append(("fit on both domains together", ens.rmsf_ca(ens.idx(u))))
        with open(os.path.join(OUT, "rmsf", f"{pid}.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["chain", "residue", "resuper_default_domain", "reference_mode_domain"]
                       + [f"rmsf_{lab.replace(' ', '_')}" for lab, _ in frames])
            for i, k in enumerate(ens.keys):
                w.writerow([k[0], k[1], per_a[i], per_b[i]] + [round(float(r[i]), 3) for _, r in frames])
        sel_rmsf = frames[0][1]
        dom_idx = ens.idx(bodies[0])
        typical = float(np.median(sel_rmsf[dom_idx]))
        dtext = ", ".join(f"{ens.keys[i][1]} ({sel_rmsf[i]:.2f})" for i in diff[:12])
        L.append(f"   {pid}: {len(diff)} residue(s) differ; median RMSF in domain 1 frame "
                 f"{typical:.2f} A; differing residues (RMSF, A): {dtext}"
                 + (" ..." if len(diff) > 12 else ""))
        tracks = [("RESUPER default", per_a), ("reference mode", per_b)]
        fig_rmsf(ens, frames, tracks, ss_ranges(pid), diff,
                 f"{pid}: Cα flexibility and domain assignments",
                 os.path.join(OUT, f"fig6_rmsf_{pid}.png"))
        if len(bodies) > 1 and pid in ("1VVD", "1VVE", "2D21", "2KR6"):
            names = [f"domain {j + 1} ({d.text})" for j, d in enumerate(bodies_rep.domains)]
            fig_ensemble(ens, [ens.idx(s) for s in bodies], names,
                         f"{pid}: all {ens.st.n_models} models, Cα traces "
                         f"(domains from {spec['bodies_from'].replace('_', ' ')})",
                         os.path.join(OUT, f"fig7_ensemble_{pid}.png"))

    # ---- 2KR6: most rigid 8-residue window inside the N-terminal segment
    ens = Ensemble("2KR6")
    repA, sA = domains("2KR6", A, ens)
    from resuper.geometry import rmsd_to_mean_converged
    seg = ens.idx(sA[0])
    win = []
    for s in range(len(seg) - 7):
        w8 = seg[s:s + 8]
        ensw = ens.X[:, w8].reshape(ens.X.shape[0], -1, 3)
        win.append((rmsd_to_mean_converged(ensw), ens.keys[w8[0]][1], ens.keys[w8[-1]][1]))
    win.sort()
    L.append(f"\n== 2KR6 N-terminal domain ({repA.domains[0].text}, {repA.domains[0].rmsd} A): "
             f"most rigid 8-residue windows: " + "; ".join(f"{a}..{b} {r:.2f} A" for r, a, b in win[:3])
             + f"; least rigid: {win[-1][1]}..{win[-1][2]} {win[-1][0]:.2f} A")
    L.append(f"   HELIX records: " + "; ".join(f"{c}{a}-{b}" for k, c, a, b in ss_ranges("2KR6") if k == "helix"))

    # ---- states detector
    rows = []
    L.append("\n== RESUPER --states (multiconformer detector)")
    for pid in ["1VVD", "1VVE", "2D21", "2KR6", "1CFC", "1X02", "1D1D", "1QO6"]:
        rc, line = run_states(pid)
        rows.append({"pdb_id": pid, "exit_code": rc, "result": line})
        L.append(f"   {pid}: {line or 'no states line (exit ' + str(rc) + ')'}")
    with open(os.path.join(OUT, "states.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb_id", "exit_code", "result"])
        w.writeheader()
        w.writerows(rows)

    text = "\n".join(L) + "\n"
    open(os.path.join(OUT, "summary.txt"), "w").write(text)
    log.info("summary:\n%s", text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
