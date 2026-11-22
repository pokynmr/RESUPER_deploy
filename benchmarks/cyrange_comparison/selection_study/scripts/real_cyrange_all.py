"""RESUPER vs REAL CYRANGE 2.0 on the 60 Phase 3 proteins and the 162 confirmation
proteins (decision S13). Real CYRANGE domains come from wwPDB validation reports
(wwpdb_cyrange.py; validated 108/108 exact against saved CYRANGE reports).

Same comparison functions as Phase 5 / S11 (compare, annotation_agreement, motion,
two-part check) and the same tests H1-H4, with real CYRANGE replacing reference
mode as the comparator. NOTE: H1-H4 were pre-registered against reference mode;
re-running them against real CYRANGE is a new, post hoc analysis (same rules).

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/real_cyrange_all.py
Outputs: selection_study/results/real_cyrange/
"""
import csv, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402
import bench_config as cfg  # noqa: E402
import phase5_compare as p5  # noqa: E402
from report_parser import parse_report, Report, Domain  # noqa: E402
from compare_methods import compare  # noqa: E402
from phase1_repo_cases import Ens  # noqa: E402
from wwpdb_cyrange import cyrange_domains, Rmsd  # noqa: E402
from confirm_run import mcnemar, write_csv  # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
R = os.path.join(S, "results")
OUT = os.path.join(R, "real_cyrange")
A, B, C = "resuper_default", "reference_mode", "real_cyrange"
SETS = {"original60": (os.path.join(R, "phase3", "candidates_annotated.csv"), os.path.join(R, "phase4")),
        "confirm162": (os.path.join(R, "confirm", "candidates_confirm.csv"), os.path.join(R, "confirm"))}


def seq_to_author(base, pid):
    raw = os.path.join(base, "annotations", "raw", f"{pid}_A.json")
    if not os.path.exists(raw):
        return {}
    m = (json.load(open(raw)).get("rcsb_polymer_entity_instance_container_identifiers") or {}).get(
        "auth_to_entity_poly_seq_mapping") or []
    out = {}
    for i, a in enumerate(m, start=1):
        try: out[i] = int(a)
        except (TypeError, ValueError): pass
    return out


def cyrange_report(pid, ens, calc):
    st, n, doms = cyrange_domains(pid)
    rep = Report()
    if st != "success":
        rep.fatal = st
        return rep, st, 0
    keys = set(ens.keys); dropped = 0
    for s in doms.values():
        s2 = {k for k in s if k in keys}; dropped += len(s - s2)
        if not s2:
            continue
        segs, byc = [], {}
        for ch, num in s2: byc.setdefault(ch, []).append(num)
        for ch in sorted(byc):
            v = sorted(byc[ch]); a = b = v[0]
            for x in v[1:]:
                if x == b + 1: b = x
                else: segs.append((ch, a, b)); a = b = x
            segs.append((ch, a, b))
        one = len(byc) == 1
        text = ", ".join(f"{'' if one else ch}{a}..{b}" for ch, a, b in segs)
        if one:
            segs = [("", a, b) for _, a, b in segs]      # same convention as single-chain reports
        rmsd = calc(pid, s2)
        rep.domains.append(Domain(segs, rmsd if rmsd is not None else -1.0, len(segs) - 1, len(s2), text))
    return rep, st, dropped


def main():
    log = cfg.setup_logging("real_cyrange_all")
    os.makedirs(OUT, exist_ok=True)
    calc = Rmsd()
    per, specific, status = [], [], []
    for tag, (cand_csv, base) in SETS.items():
        cands = list(csv.DictReader(open(cand_csv)))
        used = {r["pdb_id"] for r in csv.DictReader(open(os.path.join(R, "confirm", "per_protein.csv")))} \
            if tag == "confirm162" else {c["pdb_id"] for c in cands}
        p5.P4 = base
        for c in cands:
            pid, g = c["pdb_id"], c["group"]
            if pid not in used:
                continue
            ens = Ens(pid, os.path.join(S, "data", "raw_pdb", f"{pid}.pdb"))
            reps = {m: parse_report(open(os.path.join(base, "runs", pid, f"{m}.out")).read()) for m in (A, B)}
            reps[C], st, dropped = cyrange_report(pid, ens, calc)
            status.append({"set": tag, "pdb_id": pid, "wwpdb_cyrange": st, "residues_not_in_our_file": dropped})
            if st != "success":
                log.warning("%s: no real CYRANGE result (%s)", pid, st); continue
            chains = sorted({k[0] for k in ens.keys})
            sets = {m: [d.residue_set(ens.default_chain) for d in reps[m].domains] for m in (A, B, C)}
            rc, _ = compare(pid, reps[A], reps[C], A, C, chains, len(ens.keys))
            fc, _ = compare(pid, reps[B], reps[C], B, C, chains, len(ens.keys))
            src, ann = p5.annotation_domains(pid, ens.default_chain)
            r = {"set": tag, "pdb_id": pid, "group": g, "protein": c["protein"],
                 "ndom_resuper": len(sets[A]), "ndom_reference": len(sets[B]), "ndom_cyrange": len(sets[C]),
                 "agreement_resuper_vs_cyrange": rc["agreement"], "jaccard_resuper_vs_cyrange": rc["pooled_jaccard"],
                 "agreement_reference_vs_cyrange": fc["agreement"], "jaccard_reference_vs_cyrange": fc["pooled_jaccard"],
                 "coverage_resuper": rc["coverage_a"], "coverage_cyrange": rc["coverage_b"],
                 "ranges_resuper": " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in reps[A].domains),
                 "ranges_reference": " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in reps[B].domains),
                 "ranges_cyrange": " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in reps[C].domains),
                 "annotation_source": src, "n_annotated": len(ann)}
            for m, t in ((A, "resuper"), (B, "reference"), (C, "cyrange")):
                for k, v in p5.annotation_agreement(sets[m], ann).items():
                    r[f"{t}_{k}"] = v
            if g == "D" and "helix" in c["group_reason"]:
                a, b = (int(x) for x in c["group_reason"].split("helix ")[1].split(" ")[0].split("-"))
                amap = seq_to_author(base, pid)
                helix = {(ens.default_chain, amap[i]) for i in range(a, b + 1) if i in amap}
                for m, t in ((A, "resuper"), (C, "cyrange")):
                    r[f"{t}_finds_short_helix"] = bool(helix) and any(len(s & helix) >= 0.5 * len(helix) for s in sets[m])
            for m, other in ((A, C), (C, A)):
                for i, s in enumerate(sets[m]):
                    if any(len(s & t) / len(s | t) >= 0.5 for t in sets[other]):
                        continue
                    d = reps[m].domains[i]
                    rigid = 0 <= d.rmsd <= 1.0
                    bigger = [t for k, t in enumerate(sets[m]) if k != i]
                    mv = p5.motion(ens, s, max(bigger, key=len), log) if bigger else None
                    specific.append({"set": tag, "pdb_id": pid, "group": g, "reported_only_by": m, "domain": d.text,
                                     "rmsd": round(d.rmsd, 2), "internally_rigid": rigid,
                                     "moves": mv["supported_by_motion"] if mv else "n/a (only domain)",
                                     "passes_two_part_check": rigid and bool(mv) and mv["supported_by_motion"] is True})
            per.append(r)
            log.info("%s %s [%s] RES %d REF %d CYR %d | RES~CYR %s | REF~CYR %s", tag, pid, g, len(sets[A]),
                     len(sets[B]), len(sets[C]), rc["agreement"], fc["agreement"])
    write_csv(os.path.join(OUT, "per_protein.csv"), per)
    write_csv(os.path.join(OUT, "method_specific_domains.csv"), specific)
    write_csv(os.path.join(OUT, "wwpdb_status.csv"), status)

    L = ["RESUPER (default) vs REAL CYRANGE 2.0 (wwPDB validation reports). Post hoc re-run of the S11 tests "
         "with real CYRANGE as comparator; same rules.",
         f"Proteins with a real CYRANGE result: {len(per)} of {len(status)} "
         f"(missing: {', '.join(s['pdb_id'] + ' (' + s['wwpdb_cyrange'].strip() + ')' for s in status if s['wwpdb_cyrange'] != 'success') or 'none'})"]
    for scope, rows in (("confirm162", [r for r in per if r["set"] == "confirm162"]),
                        ("original60", [r for r in per if r["set"] == "original60"]),
                        ("all", per)):
        L.append(f"\n== {scope} (n = {len(rows)})")
        for who, key in (("RESUPER", "agreement_resuper_vs_cyrange"), ("reference mode", "agreement_reference_vs_cyrange")):
            cnt = {k: sum(r[key] == k for r in rows) for k in ("identical", "minor boundary differences", "substantial difference")}
            L.append(f"{who} vs real CYRANGE: identical {cnt['identical']}, minor {cnt['minor boundary differences']}, "
                     f"substantial {cnt['substantial difference']}; same domain count "
                     f"{sum(r['ndom_' + ('resuper' if who == 'RESUPER' else 'reference')] == r['ndom_cyrange'] for r in rows)}")
        ab = [r for r in rows if r["group"] in "AB" and r["n_annotated"]]
        for name, key, sign in (("H1 recovered annotated domains (RESUPER - CYRANGE)", "recovered", 1),
                                ("H2 missed annotated domains (CYRANGE - RESUPER)", "missed", -1)):
            d = np.array([sign * (r.get(f"resuper_{key}", 0) - r.get(f"cyrange_{key}", 0)) for r in ab])
            nz = d[d != 0]
            p = stats.wilcoxon(nz).pvalue if nz.size else float("nan")
            L.append(f"{name}: n = {len(d)}; RESUPER better {int((d > 0).sum())}, equal {int((d == 0).sum())}, "
                     f"CYRANGE better {int((d < 0).sum())}; totals RESUPER {sum(r.get('resuper_' + key, 0) for r in ab)}, "
                     f"CYRANGE {sum(r.get('cyrange_' + key, 0) for r in ab)}; Wilcoxon p = {p:.3g}")
        ce = [r for r in rows if r["group"] in "CE"]
        ao = sum(1 for r in ce if r["ndom_resuper"] >= 2 and r["ndom_cyrange"] < 2)
        bo = sum(1 for r in ce if r["ndom_cyrange"] >= 2 and r["ndom_resuper"] < 2)
        bt = sum(1 for r in ce if r["ndom_resuper"] >= 2 and r["ndom_cyrange"] >= 2)
        p, _ = mcnemar(ao, bo)
        L.append(f"H3 over-splitting in controls C+E (n = {len(ce)}): RESUPER only {ao}, CYRANGE only {bo}, both {bt}; exact McNemar p = {p:.3g}")
        dd = [r for r in rows if "resuper_finds_short_helix" in r]
        ao = sum(1 for r in dd if r["resuper_finds_short_helix"] and not r["cyrange_finds_short_helix"])
        bo = sum(1 for r in dd if r["cyrange_finds_short_helix"] and not r["resuper_finds_short_helix"])
        bt = sum(1 for r in dd if r["resuper_finds_short_helix"] and r["cyrange_finds_short_helix"])
        p, _ = mcnemar(ao, bo)
        L.append(f"H4 short helix reported (group D, n = {len(dd)}): RESUPER only {ao}, CYRANGE only {bo}, both {bt}, "
                 f"neither {len(dd) - ao - bo - bt}; exact McNemar p = {p:.3g}")
        for m in (A, C):
            sp = [x for x in specific if x["reported_only_by"] == m and (scope == "all" or x["set"] == scope)]
            L.append(f"Domains reported only by {m}: {len(sp)}; internally rigid {sum(x['internally_rigid'] for x in sp)}; "
                     f"pass two-part check {sum(x['passes_two_part_check'] for x in sp)}")
        for g in "ABDCE":
            rs = [r for r in rows if r["group"] == g]
            if rs:
                L.append(f"Group {g}: {len(rs)}; RESUPER identical to CYRANGE {sum(r['agreement_resuper_vs_cyrange'] == 'identical' for r in rs)}, "
                         f"reference mode identical to CYRANGE {sum(r['agreement_reference_vs_cyrange'] == 'identical' for r in rs)}; "
                         f"domains RES/REF/CYR {sum(r['ndom_resuper'] for r in rs)}/{sum(r['ndom_reference'] for r in rs)}/{sum(r['ndom_cyrange'] for r in rs)}")
    text = "\n".join(L) + "\n"
    open(os.path.join(OUT, "tests.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
