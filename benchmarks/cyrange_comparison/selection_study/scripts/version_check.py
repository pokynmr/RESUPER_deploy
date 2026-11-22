"""Re-run THIS repository's RESUPER with default parameters on every structure
still on disk (222 random-sample + 863 differing expanded-scan entries) and compare its domain report with the
earlier runs of the old repository. Read-only use of both repositories.
usage (from benchmarks/cyrange_comparison): python selection_study/scripts/v2_rerun_check.py"""
import os, sys, csv, subprocess, tempfile, shutil, re
from concurrent.futures import ProcessPoolExecutor
S = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
import bench_config as cfg  # noqa: E402
NEW_SRC = cfg.REPO_SRC  # the repository this folder belongs to
OUT = os.path.join(S, "results", "version_check"); os.makedirs(os.path.join(OUT, "runs"), exist_ok=True)

def old_out(pid):
    for f in (os.path.join(S, "results", "expanded", "runs", f"{pid}.out"),
              os.path.join(S, "results", "phase4", "runs", pid, "resuper_default.out"),
              os.path.join(S, "results", "confirm", "runs", pid, "resuper_default.out")):
        if os.path.exists(f): return f
def domains(text):  # domain/range lines only; state-detection lines are excluded
    return [l.strip() for l in text.splitlines() if re.search(r"range|domain", l, re.I) and not re.search(r"state|conformer", l, re.I)]
def run(item):
    pid, pdb = item
    o = os.path.join(OUT, "runs", f"{pid}.out")
    if not os.path.exists(o):
        with tempfile.TemporaryDirectory() as td:
            shutil.copyfile(pdb, os.path.join(td, f"{pid}.pdb"))
            env = dict(os.environ, PYTHONPATH=NEW_SRC, PYTHONDONTWRITEBYTECODE="1")
            p = subprocess.run([sys.executable, "-m", "resuper", f"{pid}.pdb"], cwd=td, env=env, capture_output=True, text=True, timeout=1800)
        open(o, "w").write(p.stdout if p.returncode == 0 else "ERROR\n" + p.stderr)
    new = open(o).read(); f = old_out(pid)
    if f is None: return pid, "no old run", "", ""
    old = open(f).read()
    if new.startswith("ERROR"): return pid, "new run failed", "", ""
    dn, do = domains(new), domains(old)
    return pid, "identical domains" if dn == do else "DIFFERENT", " / ".join(do), " / ".join(dn)

if __name__ == "__main__":
    items = {}
    for d in ("raw_pdb", "expanded_pdb"):
        if not os.path.isdir(os.path.join(S, "data", d)):
            continue  # structures are downloaded by phase4_screen_and_run.py / confirm_run.py / expanded_scan.py
        for f in os.listdir(os.path.join(S, "data", d)):
            if f.endswith(".pdb"): items.setdefault(f[:-4], os.path.join(S, "data", d, f))
    if not items:
        raise SystemExit("No structures on disk: run the pipeline steps that download them first (see README).")
    with ProcessPoolExecutor(8) as ex: rows = list(ex.map(run, sorted(items.items()), chunksize=4))
    with open(os.path.join(OUT, "v2_vs_v1.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["pdb_id", "status", "old_domains", "new_domains"]); w.writerows(rows)
    from collections import Counter
    print(len(rows), Counter(r[1] for r in rows))
    for r in rows:
        if r[1] != "identical domains": print(r)
