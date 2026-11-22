"""Step 2a: obtain the PDB files listed in the dataset manifest.

For every PDB ID in data/dataset_manifest_DRAFT.csv this script:
  1. downloads the deposited coordinates from RCSB (PDB format) into
     data/raw_pdb/<ID>.pdb, unless the file is already there;
  2. downloads the RCSB entry metadata (title, method, release and revision
     dates, model count) into data/raw_pdb/metadata/<ID>.json;
  3. checks that the file looks like a coordinate file (has ATOM records);
  4. records a SHA-256 checksum, so later runs can prove the input is unchanged.

Files are never overwritten. To re-download one, delete it first.

If RCSB cannot be reached, download the files yourself (e.g. from
https://www.rcsb.org/structure/1CPZ -> Download Files -> PDB Format) and pass the
folder with --import-dir. The files are copied in and checksummed the same way.
The PDB files must be named <ID>.pdb (any letter case).

usage (from benchmarks/cyrange_comparison):
  python scripts/fetch_structures.py                 # download everything missing
  python scripts/fetch_structures.py 1CPZ 2KR6       # only these IDs
  python scripts/fetch_structures.py --import-dir ~/Downloads/pdbs
Output: data/raw_pdb/, data/raw_pdb/fetch_status.csv, logs/fetch_<time>.log
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

import bench_config as cfg

MAX_TRIES = 4  # network retries with waits of 2, 4, 8 s


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def http_get(url, log):
    """Return the response body, retrying on network errors. A 404 is not retried."""
    req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT})
    for attempt in range(1, MAX_TRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise                      # entry does not exist / is obsolete
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        if attempt < MAX_TRIES:
            wait = 2 ** attempt
            log.warning("  %s: %s; retry %d in %d s", url, err, attempt, wait)
            time.sleep(wait)
    raise err


def looks_like_pdb(path):
    """True if the file contains coordinate records (guards against HTML error pages)."""
    with open(path, errors="replace") as fh:
        return any(line.startswith(("ATOM  ", "HETATM")) for line in fh)


def read_ids():
    with open(cfg.MANIFEST_DRAFT) as fh:
        return [row["pdb_id"].upper() for row in csv.DictReader(fh)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ids", nargs="*", help="PDB IDs (default: all in the manifest)")
    ap.add_argument("--import-dir", help="copy <ID>.pdb files from this folder "
                    "instead of downloading")
    ap.add_argument("--no-metadata", action="store_true",
                    help="skip the RCSB metadata download")
    args = ap.parse_args()

    log = cfg.setup_logging("fetch")
    os.makedirs(cfg.RAW_PDB_DIR, exist_ok=True)
    os.makedirs(cfg.META_DIR, exist_ok=True)
    manifest_ids = read_ids()
    ids = [i.upper() for i in args.ids] or manifest_ids
    unknown = [i for i in ids if i not in manifest_ids]
    if unknown:
        log.warning("not in the manifest (fetched anyway): %s", " ".join(unknown))

    status_rows = []
    for pdb_id in ids:
        dest = os.path.join(cfg.RAW_PDB_DIR, f"{pdb_id}.pdb")
        row = {"pdb_id": pdb_id, "source": "", "status": "", "sha256": "",
               "retrieved_utc": "", "metadata": "", "note": ""}
        try:
            if os.path.exists(dest):
                row.update(source="already present", status="ok")
            elif args.import_dir:
                cands = [f for f in os.listdir(args.import_dir)
                         if f.lower() == f"{pdb_id.lower()}.pdb"]
                if not cands:
                    raise FileNotFoundError(f"{pdb_id}.pdb not in {args.import_dir}")
                shutil.copyfile(os.path.join(args.import_dir, cands[0]), dest)
                row.update(source=f"imported from {args.import_dir}", status="ok",
                           retrieved_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"))
            else:
                url = cfg.RCSB_FILE_URL.format(pdb_id=pdb_id)
                data = http_get(url, log)
                tmp = dest + ".part"
                with open(tmp, "wb") as fh:
                    fh.write(data)
                os.replace(tmp, dest)
                row.update(source=url, status="ok",
                           retrieved_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"))
            if not looks_like_pdb(dest):
                os.replace(dest, dest + ".rejected")
                raise ValueError("file has no ATOM/HETATM records (renamed to .rejected)")
            row["sha256"] = sha256(dest)
            log.info("%s: %s (%s)", pdb_id, row["status"], row["source"])
        except urllib.error.HTTPError as e:
            row.update(status="failed", note=f"HTTP {e.code}: entry missing, obsolete, "
                       "or ID misread; check on rcsb.org")
            log.error("%s: %s", pdb_id, row["note"])
        except Exception as e:  # keep going: one bad entry must not stop the run
            row.update(status="failed", note=f"{type(e).__name__}: {e}")
            log.error("%s: %s", pdb_id, row["note"])

        if not args.no_metadata and row["status"] == "ok":
            meta = os.path.join(cfg.META_DIR, f"{pdb_id}.json")
            try:
                if not os.path.exists(meta):
                    data = http_get(cfg.RCSB_ENTRY_URL.format(pdb_id=pdb_id), log)
                    json.loads(data)                       # must be valid JSON
                    with open(meta, "wb") as fh:
                        fh.write(data)
                row["metadata"] = "ok"
            except Exception as e:
                row["metadata"] = f"failed: {type(e).__name__}: {e}"
                log.warning("%s: metadata %s", pdb_id, row["metadata"])
        status_rows.append(row)

    out = os.path.join(cfg.RAW_PDB_DIR, "fetch_status.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(status_rows[0]))
        w.writeheader()
        w.writerows(status_rows)
    n_ok = sum(r["status"] == "ok" for r in status_rows)
    log.info("done: %d ok, %d failed; status table %s", n_ok,
             len(status_rows) - n_ok, out)
    return 0 if n_ok == len(status_rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
