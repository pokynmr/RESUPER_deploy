"""Fetch candidate PDBs, run the oracle on each, and report those with multiple
'Optimal range' lines / multi-segment ranges (candidate multi-domain proteins).

usage: /usr/bin/python3 scripts/rcsb_fetch_oracle.py <id...>
"""
import os, re, subprocess, sys, urllib.request

ORACLE = os.environ.get("RESUPER_CYRANGE_HOST", "")  # ssh host with a licensed CYRANGE binary
REMOTE = "/tmp/resuper_md"
CAPDIR = "testdata"


def fetch(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "resuper"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read().decode()
    return data


def oracle(pdb_path):
    base = os.path.basename(pdb_path)
    subprocess.run(["ssh", ORACLE, f"mkdir -p {REMOTE}"], check=False)
    subprocess.run(["scp", "-q", pdb_path, f"{ORACLE}:{REMOTE}/{base}"], check=True)
    p = subprocess.run(["ssh", ORACLE, f"cd {REMOTE} && /usr/local/bin/cyrange {base}"],
                       capture_output=True, text=True)
    return p.stdout


if __name__ == "__main__":
    ids = sys.argv[1:]
    for pid in ids:
        try:
            data = fetch(pid)
        except Exception as e:
            print(f"{pid}: FETCH-ERR {e}")
            continue
        nmod = data.count("MODEL") 
        path = f"/tmp/{pid}.pdb"
        open(path, "w").write(data)
        out = oracle(path)
        lines = [l for l in out.splitlines() if "Optimal range" in l]
        print(f"\n== {pid}: models={nmod} optimal_lines={len(lines)} ==")
        for l in lines:
            print("   ", l.strip())
