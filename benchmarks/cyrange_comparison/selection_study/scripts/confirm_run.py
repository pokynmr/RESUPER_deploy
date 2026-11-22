"""Confirmation sample: screen, run, compare and test (decision S11).

Uses the Phase 4 and Phase 5 functions unchanged (screening rules S7; comparison
rules S8 with S8a/S8b), on the 166 frozen confirmation candidates, then runs the
pre-registered tests H1-H4 and the descriptive two-part rigid-domain check.

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/confirm_run.py            # everything
Outputs: selection_study/results/confirm/
  runs/<ID>/, screening.csv, annotations/, per_protein.csv, motion_tests.csv,
  method_specific_domains.csv, tests.txt
"""
import csv
import hashlib
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True

import numpy as np                                  # noqa: E402
from scipy import stats                             # noqa: E402
import bench_config as cfg                          # noqa: E402
import phase4_screen_and_run as p4                  # noqa: E402
import phase5_compare as p5                         # noqa: E402
from report_parser import parse_report              # noqa: E402
from compare_methods import compare                 # noqa: E402
from phase1_repo_cases import Ens                   # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
OUT = os.path.join(S, "results", "confirm")
A, B = "resuper_default", "reference_mode"


def frozen():
    path = os.path.join(OUT, "candidates_confirm.csv")
    d = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if d not in open(os.path.join(OUT, "FROZEN_confirm.txt")).read():
        raise SystemExit("confirmation list does not match its frozen checksum: stop")
    return list(csv.DictReader(open(path))), d


def seq_to_author(pid):
    import json
    raw = os.path.join(OUT, "annotations", "raw", f"{pid}_A.json")
    if not os.path.exists(raw):
        return {}
    m = (json.load(open(raw)).get("rcsb_polymer_entity_instance_container_identifiers") or {}).get(
        "auth_to_entity_poly_seq_mapping") or []
    out = {}
    for i, a in enumerate(m, start=1):
        try:
            out[i] = int(a)
        except (TypeError, ValueError):
            pass
    return out


def write_csv(path, rows):
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or ["note"], restval="")
        w.writeheader()
        w.writerows(rows or [{"note": "none"}])


def mcnemar(a_only, b_only):
    n = a_only + b_only
    return (stats.binomtest(a_only, n, 0.5).pvalue if n else float("nan")), n


def main():
    log = cfg.setup_logging("confirm")
    cands, digest = frozen()
    log.info("confirmation list verified (%s...), %d candidates", digest[:16], len(cands))
    resuper = cfg.import_resuper()

    # ---- screening (S7)
    screening = []
    for c in cands:
        pid = c["pdb_id"]
        path, note = p4.download(pid, log)
        row = {"pdb_id": pid, "group": c["group"]}
        if path is None:
            row.update({"class": "SOFTWARE COMPATIBILITY", "reasons": note, "severe": True})
        else:
            ins = p4.inspect(pid, "rcsb_download", path, resuper.pdb.parse_pdb, log)
            frac, names, labels = p4.in_chain_nonstandard(path)
            nch = len({r.chain for r in resuper.pdb.parse_pdb(path).residues})
            cls, why = p4.classify(ins, frac, names, labels, nch)
            row.update({"class": cls, "reasons": why, "n_models": ins.get("n_models_resuper"),
                        "n_residues": ins.get("n_residues_resuper"),
                        "severe": len(labels) >= 5 or cls == "SCIENTIFICALLY UNSUITABLE"})
        screening.append(row)
    write_csv(os.path.join(OUT, "screening.csv"), screening)
    usable = [c for c in cands if not next(s for s in screening if s["pdb_id"] == c["pdb_id"])["severe"]]
    log.info("screening: %s; severe/excluded from tests: %s",
             {k: sum(s["class"] == k for s in screening) for k in {s["class"] for s in screening}},
             [s["pdb_id"] for s in screening if s["severe"]])

    # ---- runs + repeat check
    _, fails = p4.run_all(cands, os.path.join(OUT, "runs"), log)
    with tempfile.TemporaryDirectory() as td:
        p4.run_all(cands, td, log)
        diffs = [f"{c['pdb_id']} {m}" for c in cands for m in (A, B)
                 if os.path.exists(os.path.join(OUT, "runs", c["pdb_id"], f"{m}.out")) and
                 open(os.path.join(OUT, "runs", c["pdb_id"], f"{m}.out")).read() !=
                 open(os.path.join(td, c["pdb_id"], f"{m}.out")).read()]
    log.info("run failures: %d; repeat check differences: %s", len(fails), " ".join(diffs) or "none")

    # ---- annotations (same script as Phase 4)
    subprocess.run([sys.executable, os.path.join(HERE, "domain_annotations.py"),
                    os.path.join(OUT, "annotations")] + [f"{c['pdb_id']}:A" for c in cands], check=True)
    p5.P4 = OUT                                    # annotation_domains reads <P4>/annotations/raw

    # ---- comparison (S8)
    per, motions, specific = [], [], []
    for c in usable:
        pid, g = c["pdb_id"], c["group"]
        ens = Ens(pid, os.path.join(S, "data", "raw_pdb", f"{pid}.pdb"))
        reps = {m: parse_report(open(os.path.join(OUT, "runs", pid, f"{m}.out")).read()) for m in (A, B)}
        sets = {m: [d.residue_set(ens.default_chain) for d in reps[m].domains] for m in (A, B)}
        row, pairs = compare(pid, reps[A], reps[B], A, B, sorted({k[0] for k in ens.keys}), len(ens.keys))
        src, ann = p5.annotation_domains(pid, ens.default_chain)
        r = {"pdb_id": pid, "group": g, "ndom_resuper": row["ndom_a"], "ndom_reference": row["ndom_b"],
             "agreement": row["agreement"], "pooled_jaccard": row["pooled_jaccard"],
             "coverage_resuper": row["coverage_a"], "coverage_reference": row["coverage_b"],
             "ranges_resuper": " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in reps[A].domains),
             "ranges_reference": " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in reps[B].domains),
             "annotation_source": src, "n_annotated": len(ann)}
        for m, tag in ((A, "resuper"), (B, "reference")):
            for k, v in p5.annotation_agreement(sets[m], ann).items():
                r[f"{tag}_{k}"] = v
        # H4: does each method report a domain holding >= 50% of the defining helix?
        if g == "D" and "helix" in c["group_reason"]:
            a, b = (int(x) for x in c["group_reason"].split("helix ")[1].split(" ")[0].split("-"))
            amap = seq_to_author(pid)
            helix = {(ens.default_chain, amap[i]) for i in range(a, b + 1) if i in amap}
            r["short_helix_author"] = f"{min(n for _, n in helix)}-{max(n for _, n in helix)}" if helix else ""
            for m, tag in ((A, "resuper"), (B, "reference")):
                r[f"{tag}_finds_short_helix"] = bool(helix) and any(len(s & helix) >= 0.5 * len(helix) for s in sets[m])
        # motion tests for all domain pairs of either method
        done = []
        for m in (A, B):
            for i in range(len(sets[m])):
                for j in range(i + 1, len(sets[m])):
                    key = frozenset([frozenset(sets[m][i]), frozenset(sets[m][j])])
                    if key in done:
                        continue
                    done.append(key)
                    res = p5.motion(ens, sets[m][i], sets[m][j], log)
                    if res:
                        motions.append({"pdb_id": pid, "group": g, "method": m,
                                        "domain_1": reps[m].domains[i].text, "domain_2": reps[m].domains[j].text, **res})
        # domains reported by only one method: two-part rigid-domain check
        for m, other in ((A, B), (B, A)):
            for i, s in enumerate(sets[m]):
                if any(len(s & t) / len(s | t) >= 0.5 for t in sets[other]):
                    continue
                d = reps[m].domains[i]
                rigid = 0 <= d.rmsd <= 1.0
                bigger = [t for k, t in enumerate(sets[m]) if k != i]
                mv = p5.motion(ens, s, max(bigger, key=len), log) if bigger else None
                specific.append({"pdb_id": pid, "group": g, "reported_only_by": m, "domain": d.text,
                                 "rmsd": d.rmsd, "internally_rigid": rigid,
                                 "moves": mv["supported_by_motion"] if mv else "n/a (only domain)",
                                 "passes_two_part_check": rigid and bool(mv) and mv["supported_by_motion"] is True})
        per.append(r)
        log.info("%s [%s] ndom %s/%s %s", pid, g, r["ndom_resuper"], r["ndom_reference"], r["agreement"])
    write_csv(os.path.join(OUT, "per_protein.csv"), per)
    write_csv(os.path.join(OUT, "motion_tests.csv"), motions)
    write_csv(os.path.join(OUT, "method_specific_domains.csv"), specific)

    # ---- pre-registered tests
    L = [f"Confirmation sample (S11): {len(cands)} drawn, {len(usable)} used in tests "
         f"(excluded: {', '.join(s['pdb_id'] for s in screening if s['severe']) or 'none'}). "
         f"Run failures {len(fails)}; repeat check: {'identical' if not diffs else 'DIFFERENCES ' + ' '.join(diffs)}."]
    ab = [r for r in per if r["group"] in "AB" and r["n_annotated"]]
    for name, key, sign in (("H1 recovered annotated domains (RESUPER - reference)", "recovered", 1),
                            ("H2 missed annotated domains (reference - RESUPER)", "missed", -1)):
        d = np.array([sign * (r.get(f"resuper_{key}", 0) - r.get(f"reference_{key}", 0)) for r in ab])
        nz = d[d != 0]
        p = stats.wilcoxon(nz).pvalue if nz.size else float("nan")
        L.append(f"{name}: n = {len(d)} proteins; RESUPER better {int((d > 0).sum())}, equal {int((d == 0).sum())}, "
                 f"reference better {int((d < 0).sum())}; mean difference {d.mean():+.2f}; "
                 f"totals RESUPER {sum(r.get('resuper_' + key, 0) for r in ab)}, "
                 f"reference {sum(r.get('reference_' + key, 0) for r in ab)}; Wilcoxon p = {p:.3g}")
    ce = [r for r in per if r["group"] in "CE"]
    a_only = sum(1 for r in ce if r["ndom_resuper"] >= 2 and r["ndom_reference"] < 2)
    b_only = sum(1 for r in ce if r["ndom_reference"] >= 2 and r["ndom_resuper"] < 2)
    both = sum(1 for r in ce if r["ndom_resuper"] >= 2 and r["ndom_reference"] >= 2)
    p, n = mcnemar(a_only, b_only)
    L.append(f"H3 over-splitting in controls C+E (n = {len(ce)}): RESUPER only {a_only}, reference only {b_only}, "
             f"both {both}; RESUPER total {a_only + both}, reference total {b_only + both}; exact McNemar p = {p:.3g}")
    dd = [r for r in per if r["group"] == "D" and "resuper_finds_short_helix" in r]
    a_only = sum(1 for r in dd if r["resuper_finds_short_helix"] and not r["reference_finds_short_helix"])
    b_only = sum(1 for r in dd if r["reference_finds_short_helix"] and not r["resuper_finds_short_helix"])
    both = sum(1 for r in dd if r["resuper_finds_short_helix"] and r["reference_finds_short_helix"])
    p, n = mcnemar(a_only, b_only)
    L.append(f"H4 short helix reported (group D, n = {len(dd)}): RESUPER only {a_only}, reference only {b_only}, "
             f"both {both}, neither {len(dd) - a_only - b_only - both}; exact McNemar p = {p:.3g}")
    for m in (A, B):
        sp = [x for x in specific if x["reported_only_by"] == m]
        L.append(f"Descriptive: domains reported only by {m}: {len(sp)}; internally rigid (<= 1.0 A) "
                 f"{sum(x['internally_rigid'] for x in sp)}; pass the two-part check "
                 f"{sum(x['passes_two_part_check'] for x in sp)}")
    for g in "ABDCE":
        rs = [r for r in per if r["group"] == g]
        if rs:
            L.append(f"Group {g}: {len(rs)} proteins; identical {sum(r['agreement'] == 'identical' for r in rs)}; "
                     f"domains RESUPER {sum(r['ndom_resuper'] for r in rs)}, reference {sum(r['ndom_reference'] for r in rs)}; "
                     f"mean coverage {np.mean([float(r['coverage_resuper']) for r in rs]):.3f} / "
                     f"{np.mean([float(r['coverage_reference']) for r in rs]):.3f}")
    text = "\n".join(L) + "\n"
    open(os.path.join(OUT, "tests.txt"), "w").write(text)
    log.info("results:\n%s", text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
