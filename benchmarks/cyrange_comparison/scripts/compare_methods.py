"""Phase 5: quantitative comparison of RESUPER default and reference mode.

Reads the raw reports written in Phase 4 (results/phase4_runs/<ID>/*.out) and
compares the methods protein by protein, following decisions D7 and D8
(docs/decisions_log.md).

Comparisons made
  main         resuper_default  vs  reference_mode      all 37 proteins
  fidelity     reference_mode   vs  cyrange2.0_saved    1A67, 1ADR
               resuper_default  vs  cyrange2.0_saved    1A67, 1ADR
  paper check  both methods     vs  CYRANGE RMSD values read from the paper's
               Figure 4 (Set B only; approximate, see decision D3)

The metrics, in plain terms
  Jaccard J(a,b) = |a ∩ b| / |a ∪ b|: shared residues / residues in either.
      1 = identical residue sets, 0 = no residue in common.
  Matching: each domain of one method is paired with at most one domain of
      the other so that the total Jaccard is as large as possible (Hungarian
      algorithm). This avoids pairing "domain 1 with domain 1" by position.
  Coverage = residues selected in any domain / residues in the structure.
  ΔRMSD = RMSD(method A domain) - RMSD(method B domain), matched pairs only.
  Shared-residue RMSD: the RMSD of the residues common to both domains of a
      pair, computed with RESUPER's own functions. Comparing it with each
      method's own RMSD shows whether that method's extra residues raise the
      RMSD (i.e. are more mobile than the shared core).

Outputs (results/phase5_comparison/)
  per_protein_comparison.csv   one row per protein (main comparison)
  domain_pairs.csv             one row per matched pair or unmatched domain
  cyrange_saved_comparison.csv 1A67/1ADR against the real CYRANGE 2.0 reports
  paper_fig4_crosscheck.csv    RMSDs against the paper's Figure 4 values
  rmsd_validation.csv          recomputed vs reported RMSD for every domain
  summary.txt                  aggregate numbers

usage (from benchmarks/cyrange_comparison):  python scripts/compare_methods.py
"""
import csv
import os
import statistics as st

import numpy as np
from scipy.optimize import linear_sum_assignment

import bench_config as cfg
from report_parser import parse_report

RUNS = os.path.join(cfg.RESULTS_DIR, "phase4_runs")
OUT = os.path.join(cfg.RESULTS_DIR, "phase5_comparison")
A, B = "resuper_default", "reference_mode"
NC = "not calculable"


# --------------------------------------------------------------------------
# Residue sets
# --------------------------------------------------------------------------
def domain_sets(report, chain_ids):
    """Residue sets of each domain, with (chain, number) identities.

    Single-chain reports print no chain letter; the structure's only chain is
    used instead, so that sets from different sources are comparable."""
    default_chain = chain_ids[0] if len(chain_ids) == 1 else ""
    return [d.residue_set(default_chain) for d in report.domains]


def union(sets):
    out = set()
    for s in sets:
        out |= s
    return out


def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else None


def category(j):
    if j == 1:
        return "identical"
    if j >= 0.8:
        return "strong"
    if j >= 0.5:
        return "partial"
    return "weak"


def compact(residues):
    """{(A,1),(A,2),(A,3),(A,7)} -> 'A1-3, A7' (for listing differences)."""
    if not residues:
        return ""
    out, run = [], []
    for ch, n in sorted(residues):
        if run and run[-1][0] == ch and run[-1][1] == n - 1:
            run.append((ch, n))
            continue
        if run:
            out.append(run)
        run = [(ch, n)]
    out.append(run)
    return ", ".join(f"{r[0][0]}{r[0][1]}" + (f"-{r[-1][1]}" if len(r) > 1 else "") for r in out)


# --------------------------------------------------------------------------
# RMSD with RESUPER's own functions (read-only use)
# --------------------------------------------------------------------------
class RmsdCalculator:
    def __init__(self):
        cfg.import_resuper()
        from resuper.core import _backbone_matrix
        from resuper.geometry import rmsd_to_mean_converged
        from resuper.pdb import parse_pdb
        self._bm, self._rmsd, self._parse = _backbone_matrix, rmsd_to_mean_converged, parse_pdb
        self._cache = {}

    def structure(self, pdb_id):
        if pdb_id not in self._cache:
            s = self._parse(os.path.join(cfg.RAW_PDB_DIR, f"{pdb_id}.pdb"))
            self._cache[pdb_id] = (s, {(r.chain, r.num): r.index for r in s.residues})
        return self._cache[pdb_id]

    def rmsd(self, pdb_id, residues, reference_convention=False):
        """Backbone RMSD of a residue set, as RESUPER reports it.

        reference_convention=True reproduces reference mode's atom set (the
        last residue of the structure contributes only its N atom)."""
        s, index = self.structure(pdb_id)
        idx = sorted(index[r] for r in residues if r in index)
        if len(idx) < 3:
            return None
        ens, _ = self._bm(s, s.conformers, idx, exclude_cterm=reference_convention)
        return self._rmsd(ens)


# --------------------------------------------------------------------------
# Comparison of two reports for one protein
# --------------------------------------------------------------------------
def compare(pdb_id, rep_a, rep_b, name_a, name_b, chain_ids, n_res, calc=None):
    sa, sb = domain_sets(rep_a, chain_ids), domain_sets(rep_b, chain_ids)
    pa, pb = union(sa), union(sb)
    row = {"pdb_id": pdb_id, "method_a": name_a, "method_b": name_b,
           "ndom_a": len(sa), "ndom_b": len(sb), "ndom_agree": len(sa) == len(sb)}
    # Coverage and pooled overlap
    row.update({
        "residues_a": len(pa), "residues_b": len(pb),
        "coverage_a": round(len(pa) / n_res, 4), "coverage_b": round(len(pb) / n_res, 4),
        "coverage_diff_a_minus_b": round((len(pa) - len(pb)) / n_res, 4),
        "pooled_jaccard": round(jaccard(pa, pb), 4) if jaccard(pa, pb) is not None else NC,
        "only_in_a": compact(pa - pb), "only_in_b": compact(pb - pa),
        "n_only_in_a": len(pa - pb), "n_only_in_b": len(pb - pa),
    })
    # Overlap of domains *within* one method (should be 0 for rigid domains)
    for tag, sets in (("a", sa), ("b", sb)):
        ov = set()
        for i in range(len(sets)):
            for j in range(i + 1, len(sets)):
                ov |= sets[i] & sets[j]
        row[f"within_method_overlap_{tag}"] = len(ov)

    # One-to-one matching on Jaccard (D7)
    pairs = []
    if sa and sb:
        J = np.array([[jaccard(x, y) or 0.0 for y in sb] for x in sa])
        ri, ci = linear_sum_assignment(-J)            # maximise total Jaccard
        matched = [(i, j) for i, j in zip(ri, ci) if J[i, j] > 0]
    else:
        J, matched = None, []
    ma, mb = {i for i, _ in matched}, {j for _, j in matched}
    for i, j in matched:
        da, db = rep_a.domains[i], rep_b.domains[j]
        p = {"pdb_id": pdb_id, "method_a": name_a, "method_b": name_b, "status": "matched",
             "domain_a": i + 1, "domain_b": j + 1, "ranges_a": da.text, "ranges_b": db.text,
             "n_res_a": da.n_residues, "n_res_b": db.n_residues,
             "n_shared": len(sa[i] & sb[j]), "jaccard": round(J[i, j], 4),
             "category": category(J[i, j]), "rmsd_a": da.rmsd, "rmsd_b": db.rmsd,
             "gaps_a": da.n_gaps, "gaps_b": db.n_gaps}
        if da.rmsd < 0 or db.rmsd < 0:
            p["delta_rmsd_a_minus_b"] = f"{NC}: degenerate (negative) RMSD"
        else:
            p["delta_rmsd_a_minus_b"] = round(da.rmsd - db.rmsd, 2)
        if calc is not None:
            shared = sa[i] & sb[j]
            r = calc.rmsd(pdb_id, shared)
            p["rmsd_shared_residues"] = round(r, 3) if r is not None else f"{NC}: < 3 shared residues"
            if r is not None and da.rmsd >= 0 and db.rmsd >= 0:
                p["rmsd_a_minus_shared"] = round(da.rmsd - r, 2)
                p["rmsd_b_minus_shared"] = round(db.rmsd - r, 2)
        pairs.append(p)
    for i, d in enumerate(rep_a.domains):
        if i not in ma:
            pairs.append({"pdb_id": pdb_id, "method_a": name_a, "method_b": name_b,
                          "status": f"unmatched in {name_a}", "domain_a": i + 1,
                          "ranges_a": d.text, "n_res_a": d.n_residues, "rmsd_a": d.rmsd,
                          "delta_rmsd_a_minus_b": f"{NC}: no matching domain"})
    for j, d in enumerate(rep_b.domains):
        if j not in mb:
            pairs.append({"pdb_id": pdb_id, "method_a": name_a, "method_b": name_b,
                          "status": f"unmatched in {name_b}", "domain_b": j + 1,
                          "ranges_b": d.text, "n_res_b": d.n_residues, "rmsd_b": d.rmsd,
                          "delta_rmsd_a_minus_b": f"{NC}: no matching domain"})

    # Splits / merges: a domain holding >= 50% of >= 2 domains of the other method
    events = []
    for tag, mine, other, me, them in (("a", sa, sb, name_a, name_b), ("b", sb, sa, name_b, name_a)):
        for i, x in enumerate(mine):
            held = [j + 1 for j, y in enumerate(other) if y and len(x & y) >= 0.5 * len(y)]
            if len(held) >= 2:
                events.append(f"{me} domain {i + 1} merges {them} domains {held}")
    row["split_merge"] = "; ".join(events)

    js = [p["jaccard"] for p in pairs if p["status"] == "matched"]
    row["n_matched"] = len(js)
    row["n_unmatched"] = len(pairs) - len(js)
    row["min_pair_jaccard"] = min(js) if js else NC
    row["mean_pair_jaccard"] = round(sum(js) / len(js), 4) if js else NC
    deltas = [p["delta_rmsd_a_minus_b"] for p in pairs
              if p["status"] == "matched" and isinstance(p["delta_rmsd_a_minus_b"], float)]
    row["max_abs_delta_rmsd"] = round(max(abs(d) for d in deltas), 2) if deltas else NC

    # Protein-level label (D7)
    pooled = row["pooled_jaccard"]
    if not sa and not sb:
        label = "both found no domain"
    elif row["ndom_agree"] and row["n_unmatched"] == 0 and js and min(js) == 1:
        label = "identical"
    elif (not row["ndom_agree"] or row["n_unmatched"] or not js or min(js) < 0.8
          or (pooled != NC and pooled < 0.8)):
        label = "substantial difference"
    else:
        label = "minor boundary differences"
    row["agreement"] = label
    return row, pairs


def describe(values):
    v = [x for x in values if isinstance(x, (int, float))]
    if not v:
        return "n/a"
    return (f"n={len(v)} mean {st.mean(v):+.3f} median {st.median(v):+.3f} "
            f"min {min(v):+.3f} max {max(v):+.3f}")


def write_csv(path, rows):
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(rows)


def main():
    log = cfg.setup_logging("compare")
    os.makedirs(OUT, exist_ok=True)
    manifest = {m["pdb_id"]: m for m in csv.DictReader(open(cfg.MANIFEST))}
    summary4 = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(RUNS, "run_summary.csv")))}
    calc = RmsdCalculator()

    rows, pairs, fidelity, paper, validation = [], [], [], [], []
    for pdb_id, m in manifest.items():
        pdir = os.path.join(RUNS, pdb_id)
        try:
            reps = {k: parse_report(open(os.path.join(pdir, f"{k}.out")).read()) for k in (A, B)}
        except FileNotFoundError as e:
            log.error("%s: missing Phase 4 output (%s); skipped", pdb_id, e)
            rows.append({"pdb_id": pdb_id, "agreement": f"{NC}: missing run output"})
            continue
        calc_s, _ = calc.structure(pdb_id)
        chain_ids = sorted({r.chain for r in calc_s.residues})
        n_res = len(calc_s.residues)

        # Validation (D8): recompute each reported RMSD from its own range.
        for k, rep in reps.items():
            for i, s in enumerate(domain_sets(rep, chain_ids)):
                r = calc.rmsd(pdb_id, s, reference_convention=(k == B))
                d = rep.domains[i]
                validation.append({"pdb_id": pdb_id, "method": k, "domain": i + 1,
                                   "reported_rmsd": d.rmsd,
                                   "recomputed_rmsd": round(r, 3) if r is not None else "",
                                   "agrees_to_0.01": r is not None and abs(round(r, 2) - d.rmsd) < 0.006})

        row, p = compare(pdb_id, reps[A], reps[B], A, B, chain_ids, n_res, calc)
        row.update({"paper_set": m["paper_set"], "paper_label": m["paper_label"],
                    "protein_name": m["protein_name"], "n_residues": n_res,
                    "n_models": m["n_models"], "chains": m["chains"],
                    "reference_reliability": summary4[pdb_id]["reference_reliability"],
                    "reliability_reasons": summary4[pdb_id]["reliability_reasons"],
                    "set_note": "deposited-entry stand-in (D2)" if m["paper_set"] == "A" else ""})
        for k, rep in reps.items():
            row[f"rmsd_wmean_{'a' if k == A else 'b'}"] = (
                round(sum(d.rmsd * d.n_residues for d in rep.domains) /
                      sum(d.n_residues for d in rep.domains), 3) if rep.domains else NC)
        rows.append(row)
        pairs += p
        log.info("%s: %s (ndom %s/%s, pooled J %s)", pdb_id, row["agreement"],
                 row["ndom_a"], row["ndom_b"], row["pooled_jaccard"])

        # Real CYRANGE 2.0 report (1A67, 1ADR)
        saved = os.path.join(pdir, "cyrange2.0_saved.out")
        if os.path.exists(saved):
            cy = parse_report(open(saved).read())
            for k in (A, B):
                r, _ = compare(pdb_id, reps[k], cy, k, "cyrange2.0_saved", chain_ids, n_res)
                fidelity.append(r)

        # Paper Figure 4 cross-check (Set B; approximate)
        if m["paper_set"] == "B" and m.get("paper_cyrange_rmsd_A"):
            # Paper values as a sorted list (domain order in the figure is not
            # sequence order). "Match" = same number of domains and every RMSD
            # within 0.01 A after sorting. RMSD only: the paper gives ranges as
            # drawn bars, so ranges cannot be compared.
            import re
            pv = sorted(float(x) for x in re.findall(r"\d+\.\d+", m["paper_cyrange_rmsd_A"]))

            def matches(rep):
                mv = sorted(d.rmsd for d in rep.domains)
                return len(mv) == len(pv) and all(abs(a - b) <= 0.0101 for a, b in zip(mv, pv))
            paper.append({"pdb_id": pdb_id,
                          "paper_ndom": len(pv),
                          "default_matches_paper": matches(reps[A]),
                          "reference_matches_paper": matches(reps[B]),
                          "paper_fig4_cyrange_rmsd": m["paper_cyrange_rmsd_A"],
                          "resuper_default_rmsd": " + ".join(f"{d.rmsd:.2f}" for d in reps[A].domains),
                          "reference_mode_rmsd": " + ".join(f"{d.rmsd:.2f}" for d in reps[B].domains),
                          "note": "Fig 4 values read from a figure; paper used 2010 PDB files "
                                  "and an earlier CYRANGE version"})

    ordered = ["pdb_id", "paper_set", "paper_label", "protein_name", "n_residues", "n_models",
               "chains", "agreement", "ndom_a", "ndom_b", "ndom_agree", "n_matched",
               "n_unmatched", "min_pair_jaccard", "mean_pair_jaccard", "pooled_jaccard",
               "coverage_a", "coverage_b", "coverage_diff_a_minus_b", "rmsd_wmean_a",
               "rmsd_wmean_b", "max_abs_delta_rmsd", "split_merge", "within_method_overlap_a",
               "within_method_overlap_b", "n_only_in_a", "only_in_a", "n_only_in_b", "only_in_b",
               "reference_reliability", "reliability_reasons", "set_note", "method_a", "method_b",
               "residues_a", "residues_b"]
    rows = [{k: r.get(k, "") for k in ordered} for r in rows]
    write_csv(os.path.join(OUT, "per_protein_comparison.csv"), rows)
    write_csv(os.path.join(OUT, "domain_pairs.csv"), pairs)
    write_csv(os.path.join(OUT, "cyrange_saved_comparison.csv"), fidelity)
    write_csv(os.path.join(OUT, "paper_fig4_crosscheck.csv"), paper)
    write_csv(os.path.join(OUT, "rmsd_validation.csv"), validation)

    # ---------------- summary ----------------
    L = []
    nv = sum(v["agrees_to_0.01"] for v in validation)
    L.append(f"RMSD validation: {nv}/{len(validation)} reported domain RMSDs reproduced "
             "to 0.01 A by the recomputation (required before trusting shared-residue RMSDs)")
    groups = [("all proteins", rows),
              ("Set A (deposited-entry stand-ins)", [r for r in rows if r["paper_set"] == "A"]),
              ("Set B (paper's PDB entries)", [r for r in rows if r["paper_set"] == "B"]),
              ("reliability standard", [r for r in rows if r["reference_reliability"] == "standard"]),
              ("reliability lower", [r for r in rows if r["reference_reliability"] == "lower"])]
    for name, g in groups:
        L.append(f"\n== {name}: {len(g)} proteins  (A = {A}, B = {B})")
        lab = {}
        for r in g:
            lab[r["agreement"]] = lab.get(r["agreement"], 0) + 1
        L.append("   agreement: " + ", ".join(f"{k} {v}" for k, v in sorted(lab.items())))
        L.append(f"   same domain count: {sum(r['ndom_agree'] is True for r in g)}/{len(g)}; "
                 f"domains A {sum(r['ndom_a'] for r in g)}, B {sum(r['ndom_b'] for r in g)}")
        L.append(f"   pooled Jaccard: {describe([r['pooled_jaccard'] for r in g])}")
        L.append(f"   coverage A: {describe([r['coverage_a'] for r in g])}")
        L.append(f"   coverage B: {describe([r['coverage_b'] for r in g])}")
        L.append(f"   coverage A - B: {describe([r['coverage_diff_a_minus_b'] for r in g])}; "
                 f"A higher {sum(r['coverage_diff_a_minus_b'] > 0 for r in g)}, "
                 f"equal {sum(r['coverage_diff_a_minus_b'] == 0 for r in g)}, "
                 f"B higher {sum(r['coverage_diff_a_minus_b'] < 0 for r in g)}")
        ids = {r["pdb_id"] for r in g}
        mp = [p for p in pairs if p["pdb_id"] in ids and p["status"] == "matched"]
        d = [p["delta_rmsd_a_minus_b"] for p in mp if isinstance(p["delta_rmsd_a_minus_b"], float)]
        L.append(f"   matched domain pairs: {len(mp)}; categories: " + ", ".join(
            f"{c} {sum(p['category'] == c for p in mp)}"
            for c in ("identical", "strong", "partial", "weak")))
        L.append(f"   delta RMSD (A - B) over matched pairs: {describe(d)}; "
                 f"A lower {sum(x < 0 for x in d)}, equal {sum(x == 0 for x in d)}, "
                 f"B lower {sum(x > 0 for x in d)}")
        um = [p for p in pairs if p["pdb_id"] in ids and p["status"] != "matched"]
        L.append(f"   unmatched domains: {len(um)} " + " ".join(
            f"{p['pdb_id']}({p['status'].split()[-1]})" for p in um))
    subst = [r for r in rows if r["agreement"] == "substantial difference"]
    L.append(f"\n== proteins with a substantial difference: {len(subst)}")
    for r in subst:
        L.append(f"   {r['pdb_id']} ({r['paper_label']}, reliability {r['reference_reliability']}): "
                 f"ndom {r['ndom_a']}/{r['ndom_b']}, pooled J {r['pooled_jaccard']}, "
                 f"min pair J {r['min_pair_jaccard']}, overlap within B {r['within_method_overlap_b']}"
                 + (f", {r['split_merge']}" if r["split_merge"] else ""))
    minor = [r["pdb_id"] for r in rows if r["agreement"] == "minor boundary differences"]
    L.append(f"\n== minor boundary differences: {len(minor)}: {' '.join(minor)}")
    L.append("\n== against the real CYRANGE 2.0 reports (1A67, 1ADR)")
    for f in fidelity:
        L.append(f"   {f['pdb_id']} {f['method_a']} vs CYRANGE: {f['agreement']}, "
                 f"pooled J {f['pooled_jaccard']}")
    L.append("\n== against the CYRANGE RMSDs printed in the paper's Figure 4 (Set B, "
             "approximate: read from a figure, 2010 PDB files, earlier CYRANGE)")
    for k, lbl in (("reference_matches_paper", B), ("default_matches_paper", A)):
        miss = [p["pdb_id"] for p in paper if not p[k]]
        L.append(f"   {lbl}: same domain count and all RMSDs within 0.01 A for "
                 f"{len(paper) - len(miss)}/{len(paper)}; not matching: {' '.join(miss)}")
    L.append("\nNote: 1PKT/2PNI, 1VVD/1VVE and 1CFC/1X02 are pairs of the same protein, "
             "so the 37 proteins are not 37 independent observations.")
    text = "\n".join(L) + "\n"
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text)
    log.info("summary:\n%s", text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
