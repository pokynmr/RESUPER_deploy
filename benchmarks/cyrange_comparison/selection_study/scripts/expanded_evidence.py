"""Evidence and ranking for the expanded scan (decision S14; rules fixed before the scan).

Candidates: every protein where RESUPER and REAL CYRANGE differ (domain count or
substantial difference): from the expanded scan (results/expanded/scan.csv) and from
the 222 already compared (results/real_cyrange/per_protein.csv).
For each: annotation agreement (CATH/SCOP/ECOD), domain RMSDs, two-part check of
method-specific domains, helix-split flag; then E1-E3 / X1-X3 and tiers as in S14.

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/expanded_evidence.py
Outputs: results/expanded/evidence.csv, ranked.csv, annotations/
"""
import csv, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402
import phase5_compare as p5  # noqa: E402
from report_parser import parse_report  # noqa: E402
from phase1_repo_cases import Ens, ss_ranges_file  # noqa: E402
from real_cyrange_all import cyrange_report  # noqa: E402
from wwpdb_cyrange import Rmsd  # noqa: E402
from confirm_run import write_csv  # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
R = os.path.join(S, "results")
OUT = os.path.join(R, "expanded")
RIGID = 1.0


def candidates():
    c = []
    scan = os.path.join(OUT, "scan.csv")
    for r in (csv.DictReader(open(scan)) if os.path.exists(scan) else []):
        if r.get("kept_for_evidence") == "True":
            c.append((r["pdb_id"], "expanded scan", os.path.join(S, "data", "expanded_pdb", f"{r['pdb_id']}.pdb"),
                      os.path.join(OUT, "runs", f"{r['pdb_id']}.out")))
    for r in csv.DictReader(open(os.path.join(R, "real_cyrange", "per_protein.csv"))):
        if r["ndom_resuper"] != r["ndom_cyrange"] or r["agreement_resuper_vs_cyrange"] == "substantial difference":
            base = os.path.join(R, "phase4" if r["set"] == "original60" else "confirm")
            c.append((r["pdb_id"], r["set"], os.path.join(S, "data", "raw_pdb", f"{r['pdb_id']}.pdb"),
                      os.path.join(base, "runs", r["pdb_id"], "resuper_default.out")))
    return c


def main():
    log = cfg.setup_logging("expanded_evidence")
    cands = candidates()
    log.info("candidates: %d", len(cands))
    ann_dir = os.path.join(OUT, "annotations")
    need = [f"{p}:A" for p, *_ in cands if not os.path.exists(os.path.join(ann_dir, "raw", f"{p}_A.json"))]
    for i in range(0, len(need), 200):
        subprocess.run([sys.executable, os.path.join(HERE, "domain_annotations.py"), ann_dir] + need[i:i + 200], check=False)
    p5.P4 = OUT
    calc = Rmsd()
    rows = []
    for pid, src_set, path, run in cands:
        try:
            ens = Ens(pid, path)
            rep = parse_report(open(run).read())
            crep, st, _ = cyrange_report(pid, ens, calc)
            sets = {"R": [d.residue_set(ens.default_chain) for d in rep.domains],
                    "C": [d.residue_set(ens.default_chain) for d in crep.domains]}
            src, ann = p5.annotation_domains(pid, ens.default_chain)
            aR, aC = p5.annotation_agreement(sets["R"], ann), p5.annotation_agreement(sets["C"], ann)
            rR = [d.rmsd for d in rep.domains]; rC = [d.rmsd for d in crep.domains]
            spec = {"R": [], "C": []}
            for m, o, reps in (("R", "C", rep), ("C", "R", crep)):
                for i, s in enumerate(sets[m]):
                    if any(len(s & t) / len(s | t) >= 0.5 for t in sets[o]):
                        continue
                    d = reps.domains[i]
                    rigid = 0 <= d.rmsd <= RIGID
                    bigger = [t for k, t in enumerate(sets[m]) if k != i]
                    mv = p5.motion(ens, s, max(bigger, key=len), log) if bigger else None
                    spec[m].append({"text": d.text, "rmsd": d.rmsd, "rigid": rigid,
                                    "moves": mv["supported_by_motion"] if mv else None,
                                    "pct": mv["percentile"] if mv and mv.get("percentile") != "" else 0,
                                    "two_part": rigid and bool(mv) and mv["supported_by_motion"] is True})
            # helix split: a passing RESUPER-only domain shares one HELIX record (>= 50 %) with another RESUPER domain
            hel = [set((ch, n) for n in range(a, b + 1)) for kind, ch, a, b in ss_ranges_file(path) if kind == "helix"]
            def in_helix(s): return [k for k, h in enumerate(hel) if len(s & h) >= 0.5 * len(s)]
            helix_split = False
            for x in spec["R"]:
                if not x["two_part"]:
                    continue
                xs = next(s for s, d in zip(sets["R"], rep.domains) if d.text == x["text"])
                hx = set(in_helix(xs))
                if hx and any(hx & set(in_helix(t)) for t in sets["R"] if t is not xs):
                    helix_split = True
            rec_R, rec_C = aR.get("recovered", 0), aC.get("recovered", 0)
            mix_R, mix_C = aR.get("mixed_domains", 0), aC.get("mixed_domains", 0)
            E1 = bool(ann) and (rec_R > rec_C or (rec_R == rec_C and mix_R < mix_C))
            X1 = bool(ann) and (rec_C > rec_R or (rec_R == rec_C and mix_C < mix_R) or aR.get("missed", 0) > aC.get("missed", 0))
            E2 = any(x > RIGID for x in rC) and all(0 <= x <= RIGID for x in rR)
            E3 = any(x["two_part"] for x in spec["R"])
            X2 = any(x["two_part"] for x in spec["C"])
            X3 = any(not x["rigid"] for x in spec["R"])
            ok = (E1 or E2 or E3) and not (X1 or X2 or X3)
            tier = (1 if E1 else 2 if E2 else 4 if helix_split else 3) if ok else ""
            rows.append({"pdb_id": pid, "source": src_set, "n_models": ens.st.n_models, "n_residues": len(ens.keys),
                         "ndom_resuper": len(rep.domains), "ndom_cyrange": len(crep.domains),
                         "ranges_resuper": " | ".join(f"{d.text} ({d.rmsd:.2f} Å)" for d in rep.domains),
                         "ranges_cyrange": " | ".join(f"{d.text} ({d.rmsd:.2f} Å)" for d in crep.domains),
                         "annotation_source": src, "n_annotated": len(ann),
                         "recovered_R": rec_R, "recovered_C": rec_C, "mixed_R": mix_R, "mixed_C": mix_C,
                         "resuper_only": "; ".join(f"{x['text']} ({x['rmsd']:.2f} Å, rigid {x['rigid']}, moves {x['moves']})" for x in spec["R"]),
                         "cyrange_only": "; ".join(f"{x['text']} ({x['rmsd']:.2f} Å, rigid {x['rigid']}, moves {x['moves']})" for x in spec["C"]),
                         "E1": E1, "E2": E2, "E3": E3, "X1": X1, "X2": X2, "X3": X3, "helix_split": helix_split,
                         "qualifies": ok, "tier": tier,
                         "ann_gain": (rec_R - rec_C) + 0.5 * (mix_C - mix_R),
                         "cyr_excess": round(max([x - RIGID for x in rC] + [0]), 2),
                         "motion_pct": max([x["pct"] or 0 for x in spec["R"]] + [0])})
            log.info("%s tier %s E%s X%s", pid, tier, [E1, E2, E3], [X1, X2, X3])
        except Exception as e:  # noqa: BLE001
            log.warning("%s: %s", pid, e)
            rows.append({"pdb_id": pid, "source": src_set, "qualifies": False, "error": str(e)[:200]})
    write_csv(os.path.join(OUT, "evidence.csv"), rows)
    q = [r for r in rows if r.get("qualifies")]
    q.sort(key=lambda r: (r["tier"], -r["ann_gain"], -r["cyr_excess"], -float(r["motion_pct"] or 0)))
    write_csv(os.path.join(OUT, "ranked.csv"), q)
    from collections import Counter
    log.info("qualifying: %d; by tier %s", len(q), Counter(r["tier"] for r in q))
    print(f"candidates {len(rows)}; qualifying {len(q)}; tiers {dict(Counter(r['tier'] for r in q))}")


if __name__ == "__main__":
    main()
