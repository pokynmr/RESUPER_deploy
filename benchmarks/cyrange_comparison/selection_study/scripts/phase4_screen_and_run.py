"""Selection study, Phase 4: download, screen and run both methods on the 60 candidates.

Implements decision S7 (docs/decisions_log.md), fixed before this script was run.

For every candidate in the FROZEN list (checksum verified first):
  1. download the deposited PDB-format file and RCSB metadata (never modified;
     SHA-256 recorded) -> selection_study/data/raw_pdb/
  2. inspect it (same checks as the earlier benchmark: models, chains, numbering,
     insertion codes, alternate locations, non-standard residues, backbone
     completeness, model consistency) and classify it as
     SUITABLE / SOFTWARE COMPATIBILITY problem / SCIENTIFICALLY UNSUITABLE
  3. run, exactly as a user would, with default parameters:
        resuper <ID>.pdb                      (resuper_default)
        resuper --method reference <ID>.pdb   (reference_mode)
     saving reports, warnings, exit codes and run times
  4. with --repeat-check: run everything again and confirm identical reports.

Outputs: selection_study/results/phase4/
  runs/<ID>/{resuper_default,reference_mode}.{out,stderr}, run_info.json
  screening.csv       one row per candidate: inspection results + classification
  domains.csv         one row per domain per method
  failed_or_excluded.csv
  environment.json

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/phase4_screen_and_run.py
  python selection_study/scripts/phase4_screen_and_run.py --repeat-check
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True

import bench_config as cfg                                   # noqa: E402
from fetch_structures import http_get, looks_like_pdb, sha256  # noqa: E402
from inspect_structures import inspect, scan_raw, STANDARD   # noqa: E402
from report_parser import parse_report, pooled               # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
SEL = os.path.join(S, "results", "phase3")
RAW = os.path.join(S, "data", "raw_pdb")
OUT = os.path.join(S, "results", "phase4")
METHODS = {"resuper_default": [], "reference_mode": ["--method", "reference"]}
TIMEOUT_S = 3600


def frozen_candidates():
    path = os.path.join(SEL, "candidates_selected.csv")
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if digest not in open(os.path.join(SEL, "FROZEN_selection.txt")).read():
        raise SystemExit("candidate list does not match its frozen checksum: stop")
    return list(csv.DictReader(open(path))), digest


def download(pid, log):
    """Return (path, note). Never overwrites an existing file."""
    os.makedirs(os.path.join(RAW, "metadata"), exist_ok=True)
    dest = os.path.join(RAW, f"{pid}.pdb")
    if not os.path.exists(dest):
        try:
            data = http_get(cfg.RCSB_FILE_URL.format(pdb_id=pid), log)
        except urllib.error.HTTPError as e:
            return None, f"no PDB-format file (HTTP {e.code})"
        open(dest + ".part", "wb").write(data)
        os.replace(dest + ".part", dest)
        if not looks_like_pdb(dest):
            os.replace(dest, dest + ".rejected")
            return None, "downloaded file has no coordinates"
    meta = os.path.join(RAW, "metadata", f"{pid}.json")
    if not os.path.exists(meta):
        try:
            open(meta, "wb").write(http_get(cfg.RCSB_ENTRY_URL.format(pdb_id=pid), log))
        except Exception as e:                               # metadata is optional
            log.warning("%s: metadata unavailable: %s", pid, e)
    return dest, "ok"


def in_chain_nonstandard(path):
    """Non-standard residues that are part of the polypeptide chain (model 1).

    A non-standard residue counts as in-chain only if it is peptide-bonded to a
    neighbour: its N within 2.0 Å of the previous residue's C, or its C within
    2.0 Å of the next residue's N. Ligands, sugars and side-chain adducts
    (e.g. NAG on Asn, a drug on Cys) are not in the chain and create no gap.
    Returns (fraction of chain residues, names, residue labels)."""
    import numpy as np
    res, order = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            key = (line[21], line[22:27], line[17:20].strip())
            if key not in res:
                res[key] = {"rec": line[:6], "atoms": {}}
                order.append(key)
            res[key]["atoms"][line[12:16].strip()] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])

    def bonded(a, b):     # C of a to N of b
        ca, nb = res[a]["atoms"].get("C"), res[b]["atoms"].get("N")
        return ca is not None and nb is not None and np.linalg.norm(ca - nb) < 2.0

    chain = []
    for i, k in enumerate(order):
        prev_ok = i > 0 and order[i - 1][0] == k[0] and bonded(order[i - 1], k)
        next_ok = i + 1 < len(order) and order[i + 1][0] == k[0] and bonded(k, order[i + 1])
        standard = res[k]["rec"] == "ATOM  " and k[2] in STANDARD
        if standard or prev_ok or next_ok:
            chain.append((k, standard))
    bad = [k for k, std in chain if not std]
    frac = len(bad) / len(chain) if chain else 1.0
    return frac, sorted({k[2] for k in bad}), [f"{k[2]}{k[1].strip()}" for k in bad]


def classify(row, frac, names, labels, n_chains_read):
    """Decision S7. Returns (class, reasons)."""
    unsuitable, software, caveats = [], [], []
    if frac > 0.20:
        unsuitable.append(f"{frac:.0%} of chain residues non-standard ({' '.join(names[:6])})")
    elif labels:
        severity = "SEVERE: " if len(labels) >= 5 else ""
        software.append(f"{severity}{len(labels)} in-chain non-standard residue(s) skipped by "
                        f"RESUPER ({' '.join(labels[:8])}{' ...' if len(labels) > 8 else ''}): "
                        f"each leaves an artificial gap")
    if int(row.get("n_residues_resuper") or 0) < 20:
        unsuitable.append(f"RESUPER reads only {row.get('n_residues_resuper')} residues")
    if int(row.get("models_with_different_backbone") or 0) > 0:
        unsuitable.append("models differ in backbone atoms")
    if n_chains_read > 1:
        unsuitable.append(f"{n_chains_read} chains read")
    if row.get("insertion_codes"):
        software.append("insertion codes (RESUPER merges residues sharing a number)")
    if row.get("altlocs"):
        caveats.append("alternate locations (RESUPER keeps the first)")
    if row.get("resuper_dropped"):
        caveats.append(f"RESUPER drops residues lacking N/CA/C: {row['resuper_dropped']}")
    if int(row.get("models_with_different_atoms") or 0) and not unsuitable:
        caveats.append("models differ in non-backbone atoms only")
    if unsuitable:
        return "SCIENTIFICALLY UNSUITABLE", "; ".join(unsuitable + software + caveats)
    if software:
        return "SOFTWARE COMPATIBILITY", "; ".join(software + caveats)
    return "SUITABLE", "; ".join(caveats) or "no issues"


def run_cli(args, pid, path):
    with tempfile.TemporaryDirectory() as td:
        shutil.copyfile(path, os.path.join(td, f"{pid}.pdb"))
        env = dict(os.environ, PYTHONPATH=cfg.REPO_SRC, PYTHONDONTWRITEBYTECODE="1")
        t0 = time.time()
        try:
            p = subprocess.run([sys.executable, "-m", "resuper", *args, f"{pid}.pdb"], cwd=td,
                               env=env, capture_output=True, text=True, timeout=TIMEOUT_S)
            return p.returncode, p.stdout, p.stderr, round(time.time() - t0, 1)
        except subprocess.TimeoutExpired:
            return None, "", f"timeout after {TIMEOUT_S} s", TIMEOUT_S


def run_all(cands, outdir, log):
    domains, fails = [], []
    for c in cands:
        pid = c["pdb_id"]
        path = os.path.join(RAW, f"{pid}.pdb")
        if not os.path.exists(path):
            continue
        pdir = os.path.join(outdir, pid)
        os.makedirs(pdir, exist_ok=True)
        info = {"pdb_id": pid, "input_sha256": sha256(path), "runs": {}}
        for label, args in METHODS.items():
            rc, out, err, secs = run_cli(args, pid, path)
            open(os.path.join(pdir, f"{label}.out"), "w").write(out)
            open(os.path.join(pdir, f"{label}.stderr"), "w").write(err)
            info["runs"][label] = {"command": "resuper " + " ".join(args + [f"{pid}.pdb"]),
                                   "exit_code": rc, "seconds": secs,
                                   "stdout_sha256": hashlib.sha256(out.encode()).hexdigest()}
            rep = parse_report(out)
            if rc != 0 or not rep.ok:
                fails.append({"pdb_id": pid, "group": c["group"], "stage": f"run {label}",
                              "class": "SOFTWARE COMPATIBILITY",
                              "reason": rep.fatal or err.strip()[-300:] or f"exit {rc}"})
                log.error("%s %s failed: %s", pid, label, fails[-1]["reason"])
                continue
            for i, d in enumerate(rep.domains, 1):
                domains.append({"pdb_id": pid, "group": c["group"], "method": label, "domain": i,
                                "ranges": d.text, "rmsd_A": d.rmsd, "n_gaps": d.n_gaps,
                                "n_residues": d.n_residues})
            log.info("%s %-15s %ss %d domain(s): %s", pid, label, secs, len(rep.domains),
                     " | ".join(f"{d.text} ({d.rmsd:.2f})" for d in rep.domains) or "none")
        json.dump(info, open(os.path.join(pdir, "run_info.json"), "w"), indent=1)
    return domains, fails


def write_csv(path, rows):
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or ["note"], restval="")
        w.writeheader()
        w.writerows(rows or [{"note": "none"}])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeat-check", action="store_true")
    ap.add_argument("--test", nargs="*", help="only these IDs, written to a separate test folder")
    args = ap.parse_args()
    log = cfg.setup_logging("selection_phase4")
    cands, digest = frozen_candidates()
    global OUT
    if args.test:
        cands = [c for c in cands if c["pdb_id"] in set(args.test)]
        OUT = os.path.join(OUT, "_test_run")
    log.info("frozen candidate list verified (%s...), %d candidates", digest[:16], len(cands))
    os.makedirs(OUT, exist_ok=True)

    if args.repeat_check:
        with tempfile.TemporaryDirectory() as td:
            run_all(cands, td, log)
            diffs = []
            for c in cands:
                for label in METHODS:
                    a = os.path.join(OUT, "runs", c["pdb_id"], f"{label}.out")
                    b = os.path.join(td, c["pdb_id"], f"{label}.out")
                    if os.path.exists(a) and (not os.path.exists(b) or open(a).read() != open(b).read()):
                        diffs.append(f"{c['pdb_id']} {label}")
        log.info("repeat check: differences: %s", " ".join(diffs) or "none")
        return 0 if not diffs else 1

    resuper = cfg.import_resuper()
    import numpy
    import scipy
    from resuper.core import Params
    import dataclasses
    json.dump({"python": sys.version.split()[0], "platform": platform.platform(),
               "numpy": numpy.__version__, "scipy": scipy.__version__,
               "resuper_version": resuper.__version__,
               "default_parameters": dataclasses.asdict(Params()),
               "frozen_selection_sha256": digest},
              open(os.path.join(OUT, "environment.json"), "w"), indent=1)

    screening, excluded = [], []
    for c in cands:
        pid = c["pdb_id"]
        path, note = download(pid, log)
        row = {"pdb_id": pid, "group": c["group"], "protein": c["protein"]}
        if path is None:
            row.update({"class": "SOFTWARE COMPATIBILITY", "reasons": note})
            excluded.append({"pdb_id": pid, "group": c["group"], "stage": "download",
                             "class": row["class"], "reason": note})
            screening.append(row)
            log.error("%s: %s", pid, note)
            continue
        try:
            ins = inspect(pid, "rcsb_download", path, resuper.pdb.parse_pdb, log)
        except Exception as e:
            row.update({"class": "SOFTWARE COMPATIBILITY",
                        "reasons": f"cannot be read: {type(e).__name__}: {e}"})
            excluded.append({"pdb_id": pid, "group": c["group"], "stage": "inspection",
                             "class": row["class"], "reason": row["reasons"]})
            screening.append(row)
            continue
        frac, names, labels = in_chain_nonstandard(path)
        n_chains = len({r.chain for r in resuper.pdb.parse_pdb(path).residues})
        cls, why = classify(ins, frac, names, labels, n_chains)
        row.update({k: ins.get(k, "") for k in (
            "n_models_file", "n_models_resuper", "n_residues_file_standard", "n_residues_resuper",
            "chains", "numbering", "numbering_jumps", "insertion_codes", "altlocs",
            "nonstandard_in_ATOM", "hetatm_residues", "models_with_different_atoms",
            "models_with_different_backbone", "resuper_dropped", "expdta")})
        row.update({"sequence_length_rcsb": c["length"], "in_chain_nonstandard_fraction": round(frac, 3),
                    "in_chain_nonstandard": " ".join(labels),
                    "input_sha256": sha256(path), "class": cls, "reasons": why})
        if cls == "SCIENTIFICALLY UNSUITABLE":
            excluded.append({"pdb_id": pid, "group": c["group"], "stage": "screening",
                             "class": cls, "reason": why})
        screening.append(row)
        log.info("%s [%s]: %s: %s", pid, c["group"], cls, why)

    write_csv(os.path.join(OUT, "screening.csv"), screening)
    domains, fails = run_all(cands, os.path.join(OUT, "runs"), log)
    write_csv(os.path.join(OUT, "domains.csv"), domains)
    write_csv(os.path.join(OUT, "failed_or_excluded.csv"), excluded + fails)
    counts = {}
    for r in screening:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    log.info("screening: %s; run failures: %d", counts, len(fails))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
