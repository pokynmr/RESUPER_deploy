"""Selection study, Phase 5: quantitative comparison (rules: decision S8).

A = RESUPER default, B = RESUPER reference mode (an imitation of CYRANGE 2.0,
not the CYRANGE program). For the 60 frozen candidates:

 1. A vs B with the earlier benchmark's procedure (Jaccard, one-to-one matching
    by the Hungarian algorithm, merge detection, coverage, RMSD, gaps).
 2. A and B vs the independent domain annotation (CATH/ECOD/SCOP; the source
    covering most residues): mixed / recovered / split / missed annotated domains.
 3. Relative motion of every pair of domains reported by either method, against
    a size-matched baseline (95th percentile = "supported by motion").
 4. Pattern labels (P1 RESUPER more domains, P2 reference more, P3 same count with
    different boundaries) and the evidence level of each difference.

Outputs: selection_study/results/phase5/
  per_protein.csv, domain_pairs.csv, annotation_agreement.csv, motion_tests.csv,
  summary.txt, fig_maps_<groups>.png

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/phase5_compare.py
"""
import csv
import json
import os
import statistics as stt
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True

import numpy as np                                        # noqa: E402
import bench_config as cfg                                # noqa: E402
from report_parser import parse_report                    # noqa: E402
from compare_methods import compare                       # noqa: E402
from phase1_repo_cases import Ens, size_matched_baseline  # noqa: E402
import structural_analysis as sa                          # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
P4 = os.path.join(S, "results", "phase4")
OUT = os.path.join(S, "results", "phase5")
A, B = "resuper_default", "reference_mode"
EXCLUDE_FROM_AGGREGATES = {"1WHE"}
RELATED = {"2KDO": "2L9N", "2L9N": "2KDO", "1F5Y": "6BYV", "6BYV": "1F5Y"}
CAP_LEN, CAP_GAP, MIN_BASE = 20, 10, 20   # S8a (was 40, 20)
plt = sa.plt
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


# ------------------------------------------------------------ annotations
def annotation_domains(pid, chain):
    """Domains of the classification covering most residues, as sets of
    (chain, author number). Returns (source, [sets])."""
    path = os.path.join(P4, "annotations", "raw", f"{pid}_A.json")
    if not os.path.exists(path):
        return "", []
    d = json.load(open(path))
    seqmap = (d.get("rcsb_polymer_entity_instance_container_identifiers") or {}).get(
        "auth_to_entity_poly_seq_mapping") or []
    amap = {}
    for i, a in enumerate(seqmap, start=1):
        try:
            amap[i] = int(a)
        except (TypeError, ValueError):
            pass                                       # insertion codes etc.: skip
    best, best_cov = ("", []), -1
    for src in ("CATH", "ECOD", "SCOP"):              # tie order
        doms = {}
        for n, f in enumerate(d.get("rcsb_polymer_instance_feature") or []):
            if f["type"] != src:
                continue
            key = next((v["values"][0] for v in f.get("additional_properties") or []
                        if v["name"].endswith("DOMAIN_ID") and v.get("values")), f"f{n}")
            for p in f["feature_positions"]:
                for s in range(p["beg_seq_id"], p["end_seq_id"] + 1):
                    if s in amap:
                        doms.setdefault(key, set()).add((chain, amap[s]))
        sets = [v for v in doms.values() if v]
        # S8b: a source whose own domains overlap is inconsistent; skip it
        if any(sets[i] & sets[j] for i in range(len(sets)) for j in range(i + 1, len(sets))):
            continue
        cov = len(set().union(*sets)) if sets else 0
        if cov > best_cov and sets:
            best, best_cov = (src, sets), cov
    return best


def annotation_agreement(method_sets, ann_sets):
    if not ann_sets:
        return {}
    mixed = 0
    for m in method_sets:
        parts = [len(m & a) for a in ann_sets]
        big = [x for x in parts if x >= 3 and x >= 0.10 * len(m)]
        mixed += len(big) >= 2
    rec = split = missed = 0
    for a in ann_sets:
        covers = [len(m & a) / len(a) for m in method_sets]
        best = int(np.argmax(covers)) if covers else -1
        if covers and covers[best] >= 0.5:
            m = method_sets[best]
            parts = [len(m & x) for x in ann_sets]
            if sum(1 for x in parts if x >= 3 and x >= 0.10 * len(m)) < 2:
                rec += 1
        if sum(1 for c in covers if c >= 0.2) >= 2:
            split += 1
        if not covers or max(covers) < 0.2:
            missed += 1
    return {"n_annotated": len(ann_sets), "mixed_domains": mixed, "recovered": rec,
            "split": split, "missed": missed}


# ------------------------------------------------------------ motion
_BASE = {}


def baseline(L1, gap, L2, log):
    key = (min(L1, CAP_LEN), min(max(gap, 0), CAP_GAP), min(L2, CAP_LEN))
    if key not in _BASE:
        _BASE[key] = np.array([b for _, b in size_matched_baseline(*key, log)])
    return key, _BASE[key]


def motion(ens, s1, s2, log):
    i1, i2 = ens.idx(s1), ens.idx(s2)
    if len(i1) < 3 or len(i2) < 3:
        return None
    if i1[0] > i2[0]:
        i1, i2 = i2, i1
    ang, _ = ens.rotation(i1, i2)
    gap = min(abs(a - b) for a in i1 for b in i2) - 1
    key, base = baseline(len(i1), gap, len(i2), log)
    obs = float(np.median(ang))
    if base.size < MIN_BASE:
        return {"median_deg": round(obs, 1), "baseline": f"not calculable ({base.size} pairs)",
                "percentile": "", "supported_by_motion": "not calculable"}
    pct = float((base < obs).mean() * 100)
    return {"median_deg": round(obs, 1), "baseline_key": str(key),
            "baseline_p95": round(float(np.percentile(base, 95)), 1),
            "percentile": round(pct, 1),
            "supported_by_motion": bool(obs > np.percentile(base, 95))}   # plain bool (S8a fix)


# ------------------------------------------------------------ figure
def colour_map(sets, ann):
    """One-to-one colours: each method domain is paired with at most one annotated
    domain (Hungarian matching on residue overlap) and takes its colour; unpaired
    domains get colours after the annotated ones, so no two domains of one
    method share a colour."""
    from scipy.optimize import linear_sum_assignment
    cols = {}
    if sets and ann:
        M = np.array([[len(s & a) for a in ann] for s in sets], dtype=float)
        for i, j in zip(*linear_sum_assignment(-M)):
            if M[i, j] > 0:
                cols[i] = j
    nxt = len(ann)
    for i in range(len(sets)):
        if i not in cols:
            cols[i] = nxt
            nxt += 1
    return {i: PALETTE[c % len(PALETTE)] for i, c in cols.items()}


def fig_maps(entries, path, title):
    n = len(entries)
    fig, axes = plt.subplots(n, 1, figsize=(7.0, 0.66 * n + 0.9), dpi=200)
    axes = np.atleast_1d(axes)
    for ax, e in zip(axes, entries):
        keys = e["keys"]
        x = {k: i for i, k in enumerate(keys)}
        rows = [("annotation", e["ann"]), ("RESUPER", e["A"]), ("reference", e["B"])]
        for y, (name, sets) in enumerate(rows):
            yy = 2 - y
            ax.add_patch(plt.Rectangle((0, yy - 0.32), len(keys), 0.64, color="#e2e1dc", lw=0))
            cmap = ({j: PALETTE[j % len(PALETTE)] for j in range(len(sets))} if name == "annotation"
                    else colour_map(sets, e["ann"]))
            count = {}
            for j, s in enumerate(sets):
                for k in s:
                    if k in x:
                        ax.add_patch(plt.Rectangle((x[k], yy - 0.32), 1, 0.64, color=cmap[j], lw=0))
                        count[k] = count.get(k, 0) + 1
            for k, c in count.items():            # residue in two domains of one method
                if c > 1:
                    ax.add_patch(plt.Rectangle((x[k], yy - 0.32), 1, 0.64, fill=False, hatch="////",
                                               edgecolor=sa.INK, lw=0))
            ax.text(-len(keys) * 0.01, yy, name, ha="right", va="center", fontsize=5.5, color=sa.INK2)
        ax.text(0, -0.62, str(keys[0][1]), fontsize=5, color=sa.MUTED, ha="left", va="center")
        ax.text(len(keys), -0.62, str(keys[-1][1]), fontsize=5, color=sa.MUTED, ha="right", va="center")
        ax.set_xlim(-len(keys) * 0.16, len(keys) * 1.01)
        ax.set_ylim(-0.8, 2.5)
        ax.axis("off")
        ax.set_title(f"{e['pid']} (group {e['group']}): {e['pattern']}", loc="left", fontsize=6.5,
                     color=sa.INK, pad=1.5)
    fig.suptitle(title, x=0.01, ha="left", fontsize=8, color=sa.INK, fontweight="bold")
    fig.text(0.01, 0.004, "Each colour is one domain. A method's domain is paired one-to-one with the annotated "
             "(CATH/ECOD/SCOP) domain it overlaps most and shares its colour;\nunpaired domains get new colours. "
             "Grey: in no domain. Black hatching: residue placed in two domains by one method. "
             "Numbers: first and last residue.", fontsize=5, color=sa.INK2)
    fig.tight_layout(rect=(0, 0.025, 1, 0.97))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def write_csv(path, rows):
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields and not k.startswith("_")]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    log = cfg.setup_logging("selection_phase5")
    os.makedirs(OUT, exist_ok=True)
    screening = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(P4, "screening.csv")))}
    cands = list(csv.DictReader(open(os.path.join(S, "results", "phase3", "candidates_annotated.csv"))))
    per, pairs_all, ann_rows, motions, maps = [], [], [], [], []
    for c in cands:
        pid, g = c["pdb_id"], c["group"]
        path = os.path.join(S, "data", "raw_pdb", f"{pid}.pdb")
        ens = Ens(pid, path)
        chains = sorted({k[0] for k in ens.keys})
        reps = {m: parse_report(open(os.path.join(P4, "runs", pid, f"{m}.out")).read()) for m in (A, B)}
        sets = {m: [d.residue_set(ens.default_chain) for d in reps[m].domains] for m in (A, B)}
        row, pairs = compare(pid, reps[A], reps[B], A, B, chains, len(ens.keys))
        pairs_all += pairs
        if row["ndom_a"] > row["ndom_b"]:
            pattern = "P1 RESUPER more domains"
        elif row["ndom_a"] < row["ndom_b"]:
            pattern = "P2 reference more domains"
        elif row["agreement"] == "substantial difference":
            pattern = "P3 same count, different boundaries"
        else:
            pattern = row["agreement"]
        src, ann = annotation_domains(pid, ens.default_chain)
        r = {"pdb_id": pid, "group": g, "protein": c["protein"], "n_residues": len(ens.keys),
             "n_models": ens.st.n_models, "screening": screening[pid]["class"], "pattern": pattern,
             "agreement": row["agreement"], "ndom_resuper": row["ndom_a"], "ndom_reference": row["ndom_b"],
             "pooled_jaccard": row["pooled_jaccard"], "min_pair_jaccard": row["min_pair_jaccard"],
             "coverage_resuper": row["coverage_a"], "coverage_reference": row["coverage_b"],
             "gaps_resuper": sum(d.n_gaps for d in reps[A].domains),
             "gaps_reference": sum(d.n_gaps for d in reps[B].domains),
             "rmsd_resuper": " / ".join(f"{d.rmsd:.2f}" for d in reps[A].domains),
             "rmsd_reference": " / ".join(f"{d.rmsd:.2f}" for d in reps[B].domains),
             "ranges_resuper": " | ".join(d.text for d in reps[A].domains),
             "ranges_reference": " | ".join(d.text for d in reps[B].domains),
             "split_merge": row["split_merge"], "overlap_within_reference": row["within_method_overlap_b"],
             "annotation_source": src, "n_annotated_domains": len(ann),
             "related_entry": RELATED.get(pid, ""), "literature_lead": c["literature_lead"]}
        # smaller regions found by one method only
        small_a = [p["ranges_a"] for p in pairs if p["status"] == f"unmatched in {A}" and int(p["n_res_a"]) < 40]
        small_b = [p["ranges_b"] for p in pairs if p["status"] == f"unmatched in {B}" and int(p["n_res_b"]) < 40]
        r["smaller_region_only_resuper"] = " | ".join(small_a)
        r["smaller_region_only_reference"] = " | ".join(small_b)
        for m, tag in ((A, "resuper"), (B, "reference")):
            ag = annotation_agreement(sets[m], ann)
            for k, v in ag.items():
                r[f"{tag}_{k}"] = v
            if ag:
                ann_rows.append({"pdb_id": pid, "group": g, "method": m, "source": src, **ag})
        # relative motion for every domain pair of either method
        done = []
        for m in (A, B):
            for i in range(len(sets[m])):
                for j in range(i + 1, len(sets[m])):
                    key = frozenset([frozenset(sets[m][i]), frozenset(sets[m][j])])
                    if key in done:
                        continue
                    done.append(key)
                    res = motion(ens, sets[m][i], sets[m][j], log)
                    if res is None:
                        continue
                    others = [mm for mm in (A, B) if any(set(s) == sets[m][i] for s in sets[mm])
                              and any(set(s) == sets[m][j] for s in sets[mm])]
                    motions.append({"pdb_id": pid, "group": g, "reported_by": " & ".join(others),
                                    "domain_1": reps[m].domains[i].text, "domain_2": reps[m].domains[j].text,
                                    **res})
        # evidence level for the extra split(s)
        ms = [x for x in motions if x["pdb_id"] == pid]
        only_a = [x for x in ms if x["reported_by"] == A]
        only_b = [x for x in ms if x["reported_by"] == B]
        r["resuper_only_pairs_supported_by_motion"] = f"{sum(x['supported_by_motion'] is True for x in only_a)}/{len(only_a)}"
        r["reference_only_pairs_supported_by_motion"] = f"{sum(x['supported_by_motion'] is True for x in only_b)}/{len(only_b)}"
        per.append(r)
        if pattern not in ("identical", "minor boundary differences"):
            maps.append({"pid": pid, "group": g, "pattern": pattern, "keys": ens.keys,
                         "ann": ann, "A": sets[A], "B": sets[B]})
        log.info("%s [%s] %s; annotation %s(%d); motion pairs %d", pid, g, pattern, src, len(ann), len(ms))

    write_csv(os.path.join(OUT, "per_protein.csv"), per)
    write_csv(os.path.join(OUT, "domain_pairs.csv"), pairs_all)
    write_csv(os.path.join(OUT, "annotation_agreement.csv"), ann_rows)
    write_csv(os.path.join(OUT, "motion_tests.csv"), motions)
    for grp, title in (("BA", "Groups B and A: proteins where the methods differ"),
                       ("DCE", "Groups D, C and E: proteins where the methods differ")):
        sel = [m for m in maps if m["group"] in grp]
        if sel:
            fig_maps(sel, os.path.join(OUT, f"fig_maps_{grp}.png"), title)

    # ---------------- summary
    L = ["A = RESUPER default, B = reference mode (RESUPER's imitation of CYRANGE 2.0).",
         f"Aggregates exclude {', '.join(sorted(EXCLUDE_FROM_AGGREGATES))} (severe modified-residue gaps)."]
    for g in "BADCE":
        rs = [r for r in per if r["group"] == g and r["pdb_id"] not in EXCLUDE_FROM_AGGREGATES]
        pat = {}
        for r in rs:
            pat[r["pattern"]] = pat.get(r["pattern"], 0) + 1
        L.append(f"\n== group {g} ({len(rs)} proteins): " + ", ".join(f"{k} {v}" for k, v in sorted(pat.items())))
        L.append(f"   domains: RESUPER {sum(r['ndom_resuper'] for r in rs)}, reference {sum(r['ndom_reference'] for r in rs)}; "
                 f"mean coverage RESUPER {stt.mean(float(r['coverage_resuper']) for r in rs):.3f}, "
                 f"reference {stt.mean(float(r['coverage_reference']) for r in rs):.3f}")
        withann = [r for r in rs if r["n_annotated_domains"]]
        if withann:
            for tag in ("resuper", "reference"):
                L.append(f"   vs annotation ({len(withann)} annotated, {sum(r['n_annotated_domains'] for r in withann)} domains) "
                         f"{tag:9s}: recovered {sum(r.get(tag + '_recovered', 0) for r in withann)}, "
                         f"mixed {sum(r.get(tag + '_mixed_domains', 0) for r in withann)}, "
                         f"split {sum(r.get(tag + '_split', 0) for r in withann)}, "
                         f"missed {sum(r.get(tag + '_missed', 0) for r in withann)}")
        if g in "CE":
            L.append(f"   control (should be one domain): RESUPER >= 2 domains in "
                     f"{sum(r['ndom_resuper'] >= 2 for r in rs)}, reference in {sum(r['ndom_reference'] >= 2 for r in rs)}")
        for r in rs:
            if r["pattern"] not in ("identical", "minor boundary differences"):
                L.append(f"   {r['pdb_id']}: {r['pattern']} ({r['ndom_resuper']}/{r['ndom_reference']}); "
                         f"RESUPER-only pairs moving {r['resuper_only_pairs_supported_by_motion']}, "
                         f"reference-only pairs moving {r['reference_only_pairs_supported_by_motion']}; "
                         f"annotation {r['n_annotated_domains']} dom: RESUPER rec/mix/miss "
                         f"{r.get('resuper_recovered', '-')}/{r.get('resuper_mixed_domains', '-')}/{r.get('resuper_missed', '-')}, "
                         f"reference {r.get('reference_recovered', '-')}/{r.get('reference_mixed_domains', '-')}/{r.get('reference_missed', '-')}"
                         + (f"; RESUPER-only small: {r['smaller_region_only_resuper']}" if r["smaller_region_only_resuper"] else "")
                         + (f"; reference-only small: {r['smaller_region_only_reference']}" if r["smaller_region_only_reference"] else ""))
    ok = [x for x in motions if x["supported_by_motion"] in (True, False)]
    L.append(f"\n== motion tests: {len(motions)} domain pairs; calculable {len(ok)}; "
             f"supported by motion {sum(x['supported_by_motion'] is True for x in ok)}")
    for who in (A, B, f"{A} & {B}"):
        sub = [x for x in ok if x["reported_by"] == who]
        L.append(f"   pairs reported by {who}: {len(sub)}, moving {sum(x['supported_by_motion'] is True for x in sub)}")
    text = "\n".join(L) + "\n"
    open(os.path.join(OUT, "summary.txt"), "w").write(text)
    log.info("summary:\n%s", text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
