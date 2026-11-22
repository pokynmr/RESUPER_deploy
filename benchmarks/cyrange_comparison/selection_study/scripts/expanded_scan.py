"""Expanded scan (decision S14): RESUPER (default) vs REAL CYRANGE on the whole
eligible Phase 3 pool, minus the 226 entries already drawn.

Per entry (resumable; one JSON per entry in results/expanded/rows/):
  1. real CYRANGE domains from the wwPDB validation report (none -> skipped)
  2. download PDB file (data/expanded_pdb/), S7 screening
  3. RESUPER default run (timeout 600 s), report saved in results/expanded/runs/<ID>.out
  4. compare RESUPER with CYRANGE (same compare() as before); the structure is
     deleted unless they differ (domain count or substantial difference)

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/expanded_scan.py [--workers 4] [--limit N]
"""
import argparse, csv, json, logging, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
R = os.path.join(S, "results")
OUT = os.path.join(R, "expanded")
PDB = os.path.join(S, "data", "expanded_pdb")


def work(pid):
    import phase4_screen_and_run as p4
    from report_parser import parse_report
    from compare_methods import compare
    from phase1_repo_cases import Ens
    from real_cyrange_all import cyrange_report
    from wwpdb_cyrange import cyrange_domains
    log = logging.getLogger(f"scan.{pid}")
    p4.RAW = PDB; p4.TIMEOUT_S = 600
    row = {"pdb_id": pid}
    t0 = time.time()
    try:
        st, n, doms = cyrange_domains(pid)
        row["wwpdb_cyrange"] = st.strip()
        if st != "success":
            row["status"] = "no real CYRANGE result"; return row
        path, note = p4.download(pid, log)
        if path is None:
            row["status"] = f"download: {note}"; return row
        import resuper.pdb as rp
        ins = p4.inspect(pid, "rcsb_download", path, rp.parse_pdb, log)
        frac, names, labels = p4.in_chain_nonstandard(path)
        nch = len({r.chain for r in rp.parse_pdb(path).residues})
        cls, why = p4.classify(ins, frac, names, labels, nch)
        row.update({"screening": cls, "screening_reasons": why, "n_models": ins.get("n_models_resuper"),
                    "n_residues": ins.get("n_residues_resuper"), "n_nonstandard_in_chain": len(labels)})
        if cls == "SCIENTIFICALLY UNSUITABLE" or len(labels) >= 5:
            row["status"] = "excluded (S7)"; os.remove(path); return row
        rc, out, err, secs = p4.run_cli([], pid, path)
        row["resuper_seconds"] = secs
        os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)
        open(os.path.join(OUT, "runs", f"{pid}.out"), "w").write(out)
        rep = parse_report(out)
        if rc != 0 or not rep.ok:
            row["status"] = f"RESUPER failed (rc {rc}): {(err or rep.fatal or '')[-200:]}"; os.remove(path); return row
        ens = Ens(pid, path)
        crep, _, dropped = cyrange_report(pid, ens, lambda *_: None)
        chains = sorted({k[0] for k in ens.keys})
        c, _ = compare(pid, rep, crep, "resuper_default", "real_cyrange", chains, len(ens.keys))
        row.update({"status": "ok", "ndom_resuper": c["ndom_a"], "ndom_cyrange": c["ndom_b"],
                    "agreement": c["agreement"], "pooled_jaccard": c["pooled_jaccard"],
                    "coverage_resuper": c["coverage_a"], "coverage_cyrange": c["coverage_b"],
                    "ranges_resuper": " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in rep.domains),
                    "ranges_cyrange": " | ".join(d.text for d in crep.domains),
                    "cyrange_residues_not_in_file": dropped})
        differ = c["ndom_a"] != c["ndom_b"] or c["agreement"] == "substantial difference"
        row["kept_for_evidence"] = differ
        if not differ:
            os.remove(path)
    except Exception as e:  # noqa: BLE001
        row["status"] = f"error: {type(e).__name__}: {e}"[:300]
    finally:
        row["seconds"] = round(time.time() - t0, 1)
    return row


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--limit", type=int)
    a = ap.parse_args()
    log = cfg.setup_logging("expanded_scan")
    os.makedirs(os.path.join(OUT, "rows"), exist_ok=True); os.makedirs(PDB, exist_ok=True)
    cfg.import_resuper()
    pool = [r["pdb_id"] for r in csv.DictReader(open(os.path.join(R, "phase3", "candidate_pool.csv")))]
    used = {r["pdb_id"] for r in csv.DictReader(open(os.path.join(R, "phase3", "candidates_selected.csv")))} | \
           {r["pdb_id"] for r in csv.DictReader(open(os.path.join(R, "confirm", "candidates_confirm.csv")))}
    todo = [p for p in sorted(set(pool) - used) if not os.path.exists(os.path.join(OUT, "rows", f"{p}.json"))]
    if a.limit: todo = todo[:a.limit]
    log.info("pool %d, already drawn %d, to do %d", len(set(pool)), len(used), len(todo))
    done = 0; t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(work, p): p for p in todo}
        for f in as_completed(futs):
            r = f.result(); done += 1
            json.dump(r, open(os.path.join(OUT, "rows", f"{r['pdb_id']}.json"), "w"))
            if done % 50 == 0 or done == len(todo):
                log.info("%d/%d done (%.0f s elapsed)", done, len(todo), time.time() - t0)
    rows = [json.load(open(os.path.join(OUT, "rows", f))) for f in sorted(os.listdir(os.path.join(OUT, "rows")))]
    fields = []
    for r in rows: fields += [k for k in r if k not in fields]
    with open(os.path.join(OUT, "scan.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval=""); w.writeheader(); w.writerows(rows)
    from collections import Counter
    log.info("statuses: %s", Counter(r["status"].split(":")[0] for r in rows))
    log.info("kept for evidence: %d", sum(bool(r.get("kept_for_evidence")) for r in rows))


if __name__ == "__main__":
    main()
