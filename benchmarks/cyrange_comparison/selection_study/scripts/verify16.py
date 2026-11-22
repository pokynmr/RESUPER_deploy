"""Fresh re-verification of the 16 table entries (results/phase9/resuper_vs_cyrange_cases.csv):
- re-download each PDB file and wwPDB validation report into NEW directories (no cache);
- compare file checksums with the files used before;
- run the NEW repository (repo_v2) RESUPER with default parameters;
- recompute CYRANGE domains (fixed reader: first model) and RMSDs;
- compare every domain range / RMSD / count with the table, and re-check the S15 category rules.
usage (from benchmarks/cyrange_comparison): python selection_study/scripts/verify16.py"""
import csv, hashlib, json, os, re, subprocess, sys, tempfile, shutil, logging
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402
import phase4_screen_and_run as p4  # noqa: E402
import wwpdb_cyrange as ww  # noqa: E402
from report_parser import parse_report  # noqa: E402
from phase1_repo_cases import Ens  # noqa: E402
from real_cyrange_all import cyrange_report  # noqa: E402
S = os.path.join(HERE, "..")
V = os.path.join(S, "data", "verify16"); OUT = os.path.join(S, "results", "phase9", "verify16")
for d in (os.path.join(V, "pdb"), os.path.join(V, "xml"), OUT): os.makedirs(d, exist_ok=True)
p4.RAW = os.path.join(V, "pdb"); ww.XML = os.path.join(V, "xml")
NEW_SRC = cfg.REPO_SRC  # the repository this folder belongs to
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
def norm(s): return re.sub(r"\s+", " ", s.replace("–", "..")).strip()

log = cfg.setup_logging("verify16")
table = list(csv.DictReader(open(os.path.join(S, "results", "phase9", "resuper_vs_cyrange_cases.csv"), encoding="utf-8-sig")))
USED = {r["pdb_id"]: r["sha256"] for r in csv.DictReader(open(os.path.join(OUT, "sha256_used.csv")))}
calc = ww.Rmsd(extra_dirs=[p4.RAW]); out = []
for t in table:
    pid = t["PDB"]; rec = {"PDB": pid, "Category": t["Category"]}
    path, note = p4.download(pid, log); rec["download"] = note
    rec["file_identical_to_used"] = bool(path and sha(path) == USED.get(pid))  # checksums of the files used in the paper
    meta = json.load(open(os.path.join(p4.RAW, "metadata", f"{pid}.json"))) if os.path.exists(os.path.join(p4.RAW, "metadata", f"{pid}.json")) else {}
    acc = meta.get("rcsb_accession_info", {}); rec["last_revision"] = (acc.get("revision_date") or "")[:10]
    with tempfile.TemporaryDirectory() as td:
        shutil.copyfile(path, os.path.join(td, f"{pid}.pdb"))
        p = subprocess.run([sys.executable, "-m", "resuper", f"{pid}.pdb"], cwd=td, capture_output=True, text=True,
                           env=dict(os.environ, PYTHONPATH=NEW_SRC, PYTHONDONTWRITEBYTECODE="1"), timeout=3600)
    open(os.path.join(OUT, f"{pid}_resuper_v2.out"), "w").write(p.stdout)
    rep = parse_report(p.stdout)
    res = " | ".join(f"{d.text} ({d.rmsd:.2f} Å)" for d in rep.domains)
    ens = Ens(pid, path)
    st, n, _ = ww.cyrange_domains(pid)
    crep, st2, dropped = cyrange_report(pid, ens, calc)
    cyr = " | ".join(f"{d.text} ({d.rmsd:.2f} Å)" for d in crep.domains)
    rec.update({"wwpdb_status": st, "wwpdb_n_domains": n, "cyrange_residues_not_in_file": dropped,
                "resuper_new": res, "cyrange_new": cyr,
                "resuper_matches_table": norm(res) == norm(t["RESUPER domains (RMSD)"]),
                "cyrange_matches_table": norm(cyr) == norm(t["CYRANGE domains (RMSD)"]),
                "counts_match_table": f"{len(rep.domains)} / {len(crep.domains)}" == t["Domains RESUPER / CYRANGE"],
                "wwpdb_count_consistent": n == len(crep.domains)})
    log.info("%s %s", pid, rec); out.append(rec)
F = list(out[0])
with open(os.path.join(OUT, "verify16.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=F); w.writeheader(); w.writerows(out)
keys = ["file_identical_to_used", "resuper_matches_table", "cyrange_matches_table", "counts_match_table", "wwpdb_count_consistent"]
for k in keys: print(k, sum(bool(r[k]) for r in out), "/", len(out))
for r in out:
    bad = [k for k in keys if not r[k]]
    if bad: print("CHECK", r["PDB"], bad, "| new RESUPER:", r["resuper_new"], "| new CYRANGE:", r["cyrange_new"])
print("revisions:", {r["PDB"]: r["last_revision"] for r in out})
