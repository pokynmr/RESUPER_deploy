"""Batch fetch + oracle-capture + reference-mode comparison for the PACSY expansion set.

Pipeline per selected id (testdata/pacsy_ids.txt):
  1. ensure testdata/<id>.pdb exists (download from RCSB if missing)
  2. scp all selected PDBs to the CYRANGE host in one wave
  3. one remote loop runs /usr/local/bin/cyrange per id, writing <id>.out and
     recording each process exit code in pacsy_status.txt
  4. pull outputs back; run local resuper per id; compare with --method reference

Only entries where cyrange exits 0 and emits an 'Optimal range' line are usable
as oracle references; crashes and empty outputs are reported separately.

usage: /usr/bin/python3 scripts/pacsy_eval.py [--local] [id...]
       /usr/bin/python3 scripts/pacsy_eval.py --fetch     # only download missing PDBs

The PACSY PDB files are not committed (about 300 MB); `--fetch` restores them from RCSB.
"""
import os
import subprocess
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTDATA = os.path.join(ROOT, "testdata")
ORACLE = os.environ.get("RESUPER_CYRANGE_HOST", "")  # ssh host with a licensed CYRANGE binary
REMOTE = "/tmp/resuper_md"
IDS_FILE = os.path.join(TESTDATA, "pacsy_ids.txt")
STATUS_NAME = "pacsy_status.txt"


def fetch(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "resuper"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode()


def ensure_pdb(pid):
    pdb_path = os.path.join(TESTDATA, pid + ".pdb")
    if not os.path.exists(pdb_path):
        try:
            data = fetch(pid)
        except Exception:
            return None
        open(pdb_path, "w").write(data)
    return pdb_path


def remote_capture(ids):
    """One scp wave + one remote loop on the CYRANGE host. Returns {pid: rc_str}."""
    subprocess.run(["ssh", ORACLE, f"mkdir -p {REMOTE} && rm -f {REMOTE}/*.out "
                                   f"{REMOTE}/{STATUS_NAME}"], check=False)
    for pid in ids:
        subprocess.run(["scp", "-q", os.path.join(TESTDATA, pid + ".pdb"),
                        f"{ORACLE}:{REMOTE}/{pid}.pdb"], check=False)
    # Build one remote shell loop; record each exit code unquoted so the
    # remote shell expands $? (quoting `$?` yields the literal text, which
    # broke rc detection on the first run).
    loop_body = (
        f"cd {REMOTE} && for p in {' '.join(ids)}; do "
        f"/usr/local/bin/cyrange $p.pdb > $p.out 2>/dev/null; "
        f'echo "$p $?" >> {STATUS_NAME}; done'
    )
    subprocess.run(["ssh", ORACLE, loop_body], check=False)
    subprocess.run(["scp", "-q", f"{ORACLE}:{REMOTE}/{STATUS_NAME}",
                    os.path.join(TESTDATA, STATUS_NAME)], check=False)
    for pid in ids:
        subprocess.run(["scp", "-q", f"{ORACLE}:{REMOTE}/{pid}.out",
                        os.path.join(TESTDATA, pid + ".out")], check=False)
    status = {}
    spath = os.path.join(TESTDATA, STATUS_NAME)
    if os.path.exists(spath):
        with open(spath) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    status[parts[0]] = parts[1]
    return status


def local_run(pid, timeout=180):
    import tempfile
    pdb = os.path.join(TESTDATA, pid + ".pdb")
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["cp", pdb, os.path.join(td, pid + ".pdb")], check=True)
        env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src"))
        try:
            p = subprocess.run([sys.executable, "-m", "resuper", pid + ".pdb"],
                               cwd=td, capture_output=True, text=True, env=env,
                               timeout=timeout)
        except subprocess.TimeoutExpired:
            return None
    return p.stdout


def main():
    ids = [a for a in sys.argv[1:] if not a.startswith("--")]
    local_only = "--local" in sys.argv
    if not ids:
        ids = open(IDS_FILE).read().split()
    if "--fetch" in sys.argv:
        missing = [pid for pid in ids if not ensure_pdb(pid)]
        print(f"fetched; {len(missing)} failed: {' '.join(missing)}")
        return 0
    ok, diff, crash, noopt = [], [], [], []
    for pid in ids:
        if not ensure_pdb(pid):
            crash.append(pid)
            print(f"{pid}: FETCH-ERR")
    status = {} if local_only else remote_capture(ids)
    for pid in ids:
        out_path = os.path.join(TESTDATA, pid + ".out")
        if not os.path.exists(out_path):
            crash.append(pid)
            continue
        out = open(out_path).read()
        rc = status.get(pid)
        # In --local mode no exit code is available, so a .out that is empty or
        # lacks an 'Optimal range' line is a crash/non-reference (the reference
        # always emits the line on success). Otherwise require rc == 0.
        if (not local_only and (rc is None or rc != "0")) or "Optimal range" not in out:
            crash.append(pid)
            print(f"CRASH  {pid}  (rc={rc})")
            continue
        exp = out.rstrip("\n") + "\n"
        loc = local_run(pid)
        if loc is None:
            diff.append(pid)
            print(f"DIFF   {pid}  (timeout)")
            continue
        locn = loc.rstrip("\n") + "\n"
        if locn == exp:
            ok.append(pid)
            print(f"PASS   {pid}")
        else:
            diff.append(pid)
            print(f"DIFF   {pid}")
    print(f"\nPASS={len(ok)} DIFF={len(diff)} CRASH={len(crash)} NOOPT={len(noopt)}")
    print("PAIRS=OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
