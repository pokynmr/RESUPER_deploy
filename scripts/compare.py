#!/usr/bin/env python3
"""Reference-comparison harness: run the reference `cyrange` (on $RESUPER_CYRANGE_HOST) and local
`resuper --method reference` on the same PDB + args, and compare report text.

Usage:
    compare.py testdata/final.pdb [-g 0.3] ...
    compare.py --list            # show reference cases stored in testdata/
    compare.py --capture PDB [args...]   # write oracle output to testdata/<name>.out

Both tools are run from a scratch dir containing the PDB under its base name so the
`PDB coordinate file "<name>"` header line matches.
"""
import os, re, shlex, subprocess, sys, tempfile

ORACLE = os.environ.get("RESUPER_CYRANGE_HOST", "")  # ssh host with a licensed CYRANGE binary
REMOTE_TMP = "/tmp/resuper_compare"
VERBOSE = os.environ.get("CC_VERBOSE", "0") == "1"


def norm(text: str) -> str:
    return text.rstrip("\n").replace("\r\n", "\n")


def run_local(args, cwd):
    env = dict(os.environ, PYTHONPATH=os.path.abspath("src"))
    p = subprocess.run([sys.executable, "-m", "resuper", *args],
                       cwd=cwd, capture_output=True, text=True, env=env)
    return p.returncode, p.stdout, p.stderr


def oracle_capture(args, remote: bool, cwd_name: str, local_pdb: str):
    if not remote:
        p = subprocess.run(["/usr/local/bin/cyrange", *args], cwd=cwd_name,
                           capture_output=True, text=True)
        return p.returncode, p.stdout, p.stderr
    # Push the pdb and run on the CYRANGE host, invoking with just the base filename.
    base = os.path.basename(local_pdb)
    subprocess.run(["ssh", ORACLE, f"mkdir -p {REMOTE_TMP} && rm -f {REMOTE_TMP}/*"],
                   check=False)
    subprocess.run(["scp", "-q", local_pdb, f"{ORACLE}:{REMOTE_TMP}/{base}"],
                   check=True)
    quoted = " ".join(shlex.quote(a) for a in args)
    remote_cmd = f"cd {REMOTE_TMP} && /usr/local/bin/cyrange {quoted} {shlex.quote(base)}"
    p = subprocess.run(["ssh", ORACLE, remote_cmd], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def main(argv):
    if argv and argv[0] == "--list":
        for f in sorted(os.listdir("testdata")):
            if f.endswith(".out"):
                print(f)
        return 0
    if argv and argv[0] == "--capture":
        pdb = argv[1]
        args = argv[2:]
        base = os.path.splitext(os.path.basename(pdb))[0]
        out = os.path.join("testdata", base + ".out")
        rc, stdout, _ = oracle_capture(args, True, None, pdb)
        with open(out, "w") as fh:
            fh.write(norm(stdout) + "\n")
        print(f"captured rc={rc} -> {out}")
        if VERBOSE:
            print(stdout)
        return 0
    pdb = argv[0]
    args = argv[1:]
    base = os.path.basename(pdb)
    with tempfile.TemporaryDirectory() as td:
        lpdb = os.path.join(td, base)
        subprocess.run(["cp", pdb, lpdb], check=True)
        orc = oracle_capture(args, True, td, lpdb)
        local = run_local([*args, base], td)
    if VERBOSE:
        print("=== oracle ==="); print(orc[1]); print("=== local ==="); print(local[1])
    oc, lo = norm(orc[1]), norm(local[1])
    if oc == lo:
        print(f"MATCH  ({len(oc)} bytes): {os.path.basename(pdb)} {' '.join(args)}")
        return 0
    print(f"DIFF  ({os.path.basename(pdb)} {' '.join(args)})")
    import difflib
    for line in difflib.unified_diff(oc.splitlines(), lo.splitlines(),
                                     "oracle", "local", lineterm=""):
        print(line)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
