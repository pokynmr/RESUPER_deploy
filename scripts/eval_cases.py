"""Compare local resuper output vs captured oracle .out for all main cases.

Runs the curated generalization cases, the AISAR held-out bundles
(scripts/aisar_bundles.py) and every PACSY-100 id in testdata/pacsy_ids.txt, classifies each as PASS (exact reference-mode match), DIFF
(reference captured but port disagrees) or CRASH (reference itself segfaulted,
so no usable oracle), and prints an honest aggregate on the `PASS  final` line.

usage: /usr/bin/python3 scripts/eval_cases.py [extra resuper args...]
"""
import os
import subprocess
import sys
import tempfile

EXTRA = sys.argv[1:]  # extra args, e.g. -a 2.0

# Curated generalization cases (original synthetic + single-chain solution-NMR).
CURATED = [
    "final", "m12", "m13", "m23", "m1_5", "m123", "m1234", "m10",
    "two_domain",
    "1A11", "1A6S", "1A3P", "1A67", "1d3z",
    "1a5e", "1a5r", "1a63", "1a9v", "1aca", "1adn", "1adr",
    "1aaf", "1acp", "1a90",
    "1a03", "1a7m", "1abz",
]

# Held-out real two-state ensembles from the AISAR repo (scripts/aisar_bundles.py).
AISAR = [
    "gluc_s1", "gluc_s2", "gluc_s12", "gluc_nmr", "gluc_a5", "gluc_b3", "gluc_a5b3",
    "gluc_c2", "gluc_d4", "cdk_s1", "cdk_s2", "cdk_s12", "cdk_nmr", "cdk_a1b2", "cdk_c4d3",
]


def pacsy_ids():
    path = os.path.join("testdata", "pacsy_ids.txt")
    if not os.path.exists(path):
        return []
    return [line.strip() for line in open(path) if line.strip()]


def classify(base):
    """Return ('PASS'|'DIFF'|'CRASH'|'MISSING', base)."""
    out = os.path.join("testdata", base + ".out")
    pdb = os.path.join("testdata", base + ".pdb")
    if not os.path.exists(out):
        return ("CRASH", base)
    if not os.path.exists(pdb):
        return ("MISSING", base)  # pacsy_eval.py --fetch / aisar_bundles.py --fetch --testdata
    text = open(out).read()
    # Reference segfaulted (no result line after the header) -> no usable oracle.
    if "Optimal range" not in text and "WARNING" not in text:
        return ("CRASH", base)
    expected = text.rstrip("\n") + "\n"
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["cp", pdb, os.path.join(td, os.path.basename(pdb))], check=True)
        env = dict(os.environ, PYTHONPATH=os.path.abspath("src"))
        p = subprocess.run([sys.executable, "-m", "resuper", *EXTRA,
                            os.path.basename(pdb)],
                           cwd=td, capture_output=True, text=True, env=env)
    got = p.stdout.rstrip("\n") + "\n"
    return ("PASS" if got == expected else "DIFF", base)


def main():
    # Union of curated + PACSY, dedup by id (case-sensitive); preserve the
    # curated ordering and then append any PACSY-only ids.
    seen = set()
    cases = []
    for base in CURATED + AISAR:
        if base not in seen:
            seen.add(base)
            cases.append(base)
    for base in pacsy_ids():
        if base not in seen:
            seen.add(base)
            cases.append(base)

    counts = {"PASS": 0, "DIFF": 0, "CRASH": 0, "MISSING": 0}
    for base in cases:
        status, _ = classify(base)
        counts[status] += 1
        print(f"{status}  {base}")

    total = len(cases) - counts["CRASH"] - counts["MISSING"]
    print(f"\nPASS  final  {counts['PASS']}/{total} exact match; "
          f"{counts['DIFF']} differ; {counts['CRASH']} no reference "
          f"(oracle segfault); {counts['MISSING']} PDB not fetched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
