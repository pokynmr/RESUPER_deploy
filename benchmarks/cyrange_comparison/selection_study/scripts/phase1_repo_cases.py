"""Selection study, Phase 1: the repository structures where RESUPER and CYRANGE differ.

Cases were found by screening all 110 real ensembles that have a saved CYRANGE 2.0
report in the RESUPER repository (testdata/oracle_perturb/<ID>__ident.out),
using the earlier benchmark's domain matching. Included: every ensemble where
RESUPER default and CYRANGE differ in domain count, have an unmatched or weakly
matched domain, or a split/merge, in EITHER direction:
    1A11, 11TG, 11TJ, 1DX7  (RESUPER reports more or different domains)
    6EMR                     (CYRANGE reports an extra domain; opposite direction)
    1BFW                     (degenerate CYRANGE output; RESUPER reports none)

For each case this script
  1. runs `resuper <ID>.pdb` and `resuper --method reference <ID>.pdb` afresh
     (unmodified RESUPER, default parameters) and stores both reports next to
     the saved real CYRANGE 2.0 report;
  2. compares all three pairwise (domain matching by Jaccard, as before);
  3. measures structural evidence:
     - rotation of RESUPER domain 2 relative to domain 1 in each model, judged
       against a SIZE-MATCHED baseline: segment pairs with the same lengths and
       spacing taken from inside 29 proteins that both methods call one rigid
       domain (earlier benchmark), because small segments give noisier angles;
     - for single-helix cases, the angle between the axes of the two halves in
       each model (constant = rigid curved helix; varying = hinge);
     - RMSD of CYRANGE's fused range vs RESUPER's separate domains, and the
       spread of the Cα–Cα distance between the domain centres;
     - per-residue Cα RMSF, RESUPER's states detector;
     - PDB revision history (could coordinates have changed since capture?).

Outputs: selection_study/results/phase1/
  <ID>/resuper_default.out, reference_mode.out, cyrange2.0_saved.out (+ .stderr)
  comparisons.csv, evidence.csv, baseline_size_matched.csv, summary.txt
  fig_<ID>_ensemble.png, fig_<ID>_rmsf.png, fig_rotation_vs_baseline.png

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/phase1_repo_cases.py
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True

import numpy as np                                   # noqa: E402
import bench_config as cfg                           # noqa: E402
from report_parser import parse_report               # noqa: E402
from compare_methods import compare                  # noqa: E402
import structural_analysis as sa                     # noqa: E402
from resuper.geometry import rmsd_to_mean_converged  # noqa: E402

OUT = os.path.join(cfg.BENCH_DIR, "selection_study", "results", "phase1")
CASES = ["1A11", "11TG", "11TJ", "1DX7", "6EMR", "1BFW"]
plt = sa.plt


def pdb_path(pid):
    repo = os.path.join(cfg.REPO_TESTDATA, f"{pid}.pdb")
    if os.path.exists(repo):
        return repo, "repository testdata (committed)"
    return os.path.join(cfg.DATA_DIR, "fidelity_pdb", f"{pid}.pdb"), \
        "RCSB download (earlier Phase 3); CYRANGE report in repository"


class Ens(sa.Ensemble):
    """structural_analysis.Ensemble, but read from an explicit file path."""
    def __init__(self, pid, path):                   # noqa: D401 - same fields as parent
        self.id = pid
        self.st = sa.parse_pdb(path)
        self.keys = [(r.chain, r.num) for r in self.st.residues]
        self.index = {k: i for i, k in enumerate(self.keys)}
        X = np.zeros((self.st.n_models, len(self.keys), 3, 3))
        for k, conf in enumerate(self.st.conformers):
            for i, key in enumerate(self.keys):
                for j, a in enumerate(sa.BACKBONE):
                    X[k, i, j] = conf[key][a]
        self.X = X
        chains = sorted({c for c, _ in self.keys})
        self.default_chain = chains[0] if len(chains) == 1 else ""


def ss_ranges_file(path):
    """Helix and strand residue ranges from the HELIX / SHEET records of a PDB file."""
    out = []
    for line in open(path):
        try:
            if line.startswith("HELIX "):
                out.append(("helix", line[19].strip() or "A", int(line[21:25]), int(line[33:37])))
            elif line.startswith("SHEET "):
                out.append(("strand", line[21].strip() or "A", int(line[22:26]), int(line[33:37])))
        except ValueError:
            pass
    return out


def run_cli(args, pid, path):
    with tempfile.TemporaryDirectory() as td:
        shutil.copyfile(path, os.path.join(td, f"{pid}.pdb"))
        env = dict(os.environ, PYTHONPATH=cfg.REPO_SRC, PYTHONDONTWRITEBYTECODE="1")
        p = subprocess.run([sys.executable, "-m", "resuper", *args, f"{pid}.pdb"], cwd=td,
                           env=env, capture_output=True, text=True, timeout=3600)
    return p.returncode, p.stdout, p.stderr


def helix_axis(ca):
    """Unit vector along a helix segment: first principal axis of its Cα atoms,
    oriented from the first to the last residue."""
    c = ca - ca.mean(axis=0)
    v = np.linalg.svd(c)[2][0]
    return v if np.dot(v, ca[-1] - ca[0]) > 0 else -v


def size_matched_baseline(L1, gap, L2, log):
    """Median rotation angle (per segment pair) for pairs of segments with lengths
    L1, L2 separated by `gap` residues, taken from inside single rigid domains."""
    comp = list(csv.DictReader(open(os.path.join(cfg.RESULTS_DIR, "phase5_comparison",
                                                  "per_protein_comparison.csv"))))
    meds = []
    for c in comp:
        if not (c["ndom_a"] == "1" and c["ndom_b"] == "1"):
            continue
        ens = sa.Ensemble(c["pdb_id"])
        rep = parse_report(open(os.path.join(cfg.RESULTS_DIR, "phase4_runs", c["pdb_id"],
                                             "resuper_default.out")).read())
        dom = ens.idx(rep.domains[0].residue_set(ens.default_chain))
        # contiguous stretches inside the domain
        runs, cur = [], [dom[0]]
        for i in dom[1:]:
            if i == cur[-1] + 1:
                cur.append(i)
            else:
                runs.append(cur)
                cur = [i]
        runs.append(cur)
        span = L1 + gap + L2
        for run in runs:
            for s in range(0, len(run) - span + 1, max(1, span // 2)):
                i1 = run[s:s + L1]
                i2 = run[s + L1 + gap:s + span]
                ang, _ = ens.rotation(i1, i2)
                meds.append((c["pdb_id"], float(np.median(ang))))
    log.info("size-matched baseline L1=%d gap=%d L2=%d: %d segment pairs", L1, gap, L2, len(meds))
    return meds


def main():
    log = cfg.setup_logging("selection_phase1")
    os.makedirs(OUT, exist_ok=True)
    comps, evidence, L, base_rows = [], [], [], []
    rot_plot = []
    for pid in CASES:
        path, source = pdb_path(pid)
        pdir = os.path.join(OUT, pid)
        os.makedirs(pdir, exist_ok=True)
        reps = {}
        for label, args in (("resuper_default", []), ("reference_mode", ["--method", "reference"])):
            rc, out, err = run_cli(args, pid, path)
            open(os.path.join(pdir, f"{label}.out"), "w").write(out)
            open(os.path.join(pdir, f"{label}.stderr"), "w").write(err)
            reps[label] = parse_report(out)
            if rc != 0 or not reps[label].ok:
                log.error("%s %s failed: rc=%s %s", pid, label, rc, reps[label].fatal)
        saved = os.path.join(cfg.REPO_TESTDATA, "oracle_perturb", f"{pid}__ident.out")
        shutil.copyfile(saved, os.path.join(pdir, "cyrange2.0_saved.out"))
        reps["cyrange2.0_saved"] = parse_report(open(saved).read())

        ens = Ens(pid, path)
        chains = sorted({c for c, _ in ens.keys})
        n = len(ens.keys)
        for a, b in (("resuper_default", "cyrange2.0_saved"), ("reference_mode", "cyrange2.0_saved"),
                     ("resuper_default", "reference_mode")):
            row, _ = compare(pid, reps[a], reps[b], a, b, chains, n)
            comps.append({k: row[k] for k in ("pdb_id", "method_a", "method_b", "ndom_a", "ndom_b",
                                               "pooled_jaccard", "coverage_a", "coverage_b",
                                               "split_merge", "agreement", "only_in_a", "only_in_b")})

        # ---- revision history
        meta = os.path.join(cfg.DATA_DIR, "fidelity_pdb", "metadata", f"{pid}.json")
        rev = ""
        if os.path.exists(meta):
            m = json.load(open(meta))
            h = m.get("pdbx_audit_revision_history", [])
            rev = (f"latest {m['rcsb_accession_info']['revision_date'][:10]}, highest major "
                   f"revision {max(r['major_revision'] for r in h)}" if h else "")
        ev = {"pdb_id": pid, "input": source, "n_models": ens.st.n_models, "n_residues": n,
              "chains": " ".join(chains), "revision": rev or "repository copy",
              "helix_records": " ".join(f"{a}-{b}" for k, c, a, b in ss_ranges_file(path) if k == "helix"),
              "strand_records": sum(1 for k, *_ in ss_ranges_file(path) if k == "strand")}
        for lab, rep in reps.items():
            ev[f"{lab}"] = " | ".join(f"{d.text} ({d.rmsd:.2f} A)" for d in rep.domains) or "none"

        rd = reps["resuper_default"]
        dsets = [d.residue_set(ens.default_chain) for d in rd.domains]
        if len(dsets) >= 2:
            # use the two RESUPER domains that are adjacent in sequence and largest
            i1, i2 = ens.idx(dsets[0]), ens.idx(dsets[1])
            ang, shift = ens.rotation(i1, i2)
            ev["rotation_median_deg"] = round(float(np.median(ang)), 1)
            ev["rotation_max_deg"] = round(float(ang.max()), 1)
            gap = i2[0] - i1[-1] - 1
            base = size_matched_baseline(len(i1), max(gap, 0), len(i2), log)
            bm = np.array([b for _, b in base])
            ev["baseline_n_pairs"] = len(bm)
            ev["baseline_median_deg"] = round(float(np.median(bm)), 1)
            ev["baseline_p95_deg"] = round(float(np.percentile(bm, 95)), 1)
            ev["percentile_of_observed_in_baseline"] = round(float((bm < np.median(ang)).mean() * 100), 1)
            for p_, b in base:
                base_rows.append({"case": pid, "source_protein": p_, "L1": len(i1), "gap": gap,
                                  "L2": len(i2), "median_deg": round(b, 2)})
            rot_plot.append((pid, ang, bm))
            # fused vs separate
            union = sorted(set(i1) | set(i2) | set(range(i1[-1] + 1, i2[0])))
            ensu = ens.X[:, union].reshape(ens.X.shape[0], -1, 3)
            ev["rmsd_fused_span_A"] = round(rmsd_to_mean_converged(ensu), 2)
            ev["rmsd_domains_A"] = " / ".join(f"{d.rmsd:.2f}" for d in rd.domains[:2])
            # distance between the two domain centres (Cα), per model
            ca = ens.X[:, :, 1]
            dist = np.linalg.norm(ca[:, i1].mean(axis=1) - ca[:, i2].mean(axis=1), axis=1)
            ev["centre_distance_A"] = f"{dist.min():.1f}-{dist.max():.1f} (SD {dist.std():.2f})"
            # helix bend between the two halves (only meaningful inside one helix)
            bends = [float(np.degrees(np.arccos(np.clip(np.dot(helix_axis(ca[k, i1]),
                                                               helix_axis(ca[k, i2])), -1, 1))))
                     for k in range(ca.shape[0])]
            mean_ca = ens.fit(union)[:, :, 1].mean(axis=0)
            mean_bend = float(np.degrees(np.arccos(np.clip(np.dot(helix_axis(mean_ca[i1]),
                                                                  helix_axis(mean_ca[i2])), -1, 1))))
            ev["axis_angle_per_model_deg"] = f"{min(bends):.0f}-{max(bends):.0f} (mean {np.mean(bends):.0f}, SD {np.std(bends):.1f})"
            ev["axis_angle_mean_structure_deg"] = round(mean_bend, 1)
        # states detector
        rc, out, _ = run_cli(["--states"], pid, path)
        ev["states"] = next((l.strip() for l in out.splitlines() if "Conformational states" in l), f"exit {rc}")
        evidence.append(ev)

        # figures
        rsets = [d.residue_set(ens.default_chain) for d in reps["reference_mode"].domains]
        csets = [d.residue_set(ens.default_chain) for d in reps["cyrange2.0_saved"].domains if d.rmsd >= 0]
        if dsets:
            frames = [(f"fit on RESUPER domain {j + 1}", ens.rmsf_ca(ens.idx(s))) for j, s in enumerate(dsets)]
            tracks = [("RESUPER default", sa.assignment(ens, dsets)),
                      ("reference mode", sa.assignment(ens, rsets)),
                      ("CYRANGE 2.0 (saved)", sa.assignment(ens, csets))]
            per_a, per_c = tracks[0][1], tracks[2][1]
            diff = [i for i in range(n) if (per_a[i] is None) != (per_c[i] is None)]
            sa.fig_rmsf(ens, frames, tracks, ss_ranges_file(path),
                        diff, f"{pid}: Cα flexibility; RESUPER, reference mode and real CYRANGE",
                        os.path.join(OUT, f"fig_{pid}_rmsf.png"))
            if len(dsets) >= 2:
                names = [f"RESUPER domain {j + 1} ({d.text})" for j, d in enumerate(rd.domains[:2])]
                sa.fig_ensemble(ens, [ens.idx(s) for s in dsets[:2]], names,
                                f"{pid}: all {ens.st.n_models} models, Cα traces",
                                os.path.join(OUT, f"fig_{pid}_ensemble.png"))
        log.info("%s done", pid)

    # rotation vs baseline figure
    if rot_plot:
        fig, ax = plt.subplots(figsize=(7.0, 0.55 * len(rot_plot) + 1.2), dpi=200)
        for y, (pid, ang, bm) in enumerate(rot_plot):
            lo, hi = np.percentile(bm, [5, 95])
            ax.plot([lo, hi], [y, y], color="#c9c8c2", lw=7, solid_capstyle="butt", zorder=1)
            ax.plot([np.median(bm)] * 2, [y - 0.22, y + 0.22], color=sa.INK2, lw=1, zorder=2)
            ax.scatter(ang, np.full(len(ang), y) + np.random.default_rng(y).uniform(-0.12, 0.12, len(ang)),
                       s=9, color=sa.DOMAIN[0], lw=0, alpha=0.7, zorder=3)
            ax.plot([np.median(ang)] * 2, [y - 0.3, y + 0.3], color=sa.INK, lw=1.8, zorder=4)
        ax.set_yticks(range(len(rot_plot)))
        ax.set_yticklabels([p for p, _, _ in rot_plot])
        ax.set_ylim(len(rot_plot) - 0.5, -0.8)
        ax.text(0.99, -0.6, "grey bar: 5–95% of size-matched baseline (segments inside rigid domains)",
                fontsize=6.5, color=sa.INK2, ha="right", transform=ax.get_yaxis_transform())
        ax.set_xlabel("Rotation of RESUPER domain 2 relative to domain 1, per model (degrees); black bar = median")
        for s in ("right", "top"):
            ax.spines[s].set_visible(False)
        ax.set_title("Do RESUPER's two domains move relative to each other?", loc="left",
                     fontsize=8.5, color=sa.INK, fontweight="bold")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "fig_rotation_vs_baseline.png"), facecolor="white",
                    bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)

    for name, rows in (("comparisons.csv", comps), ("evidence.csv", evidence),
                       ("baseline_size_matched.csv", base_rows)):
        fields = []
        for r in rows:
            fields += [k for k in r if k not in fields]
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, restval="")
            w.writeheader()
            w.writerows(rows)
    for ev in evidence:
        L.append(f"== {ev['pdb_id']}")
        for k, v in ev.items():
            if k != "pdb_id":
                L.append(f"   {k}: {v}")
    L.append("\n== comparisons")
    for c in comps:
        L.append(f"   {c['pdb_id']} {c['method_a']} vs {c['method_b']}: ndom {c['ndom_a']}/{c['ndom_b']}, "
                 f"pooled J {c['pooled_jaccard']}, {c['agreement']}"
                 + (f"; {c['split_merge']}" if c["split_merge"] else ""))
    text = "\n".join(L) + "\n"
    open(os.path.join(OUT, "summary.txt"), "w").write(text)
    log.info("summary:\n%s", text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
