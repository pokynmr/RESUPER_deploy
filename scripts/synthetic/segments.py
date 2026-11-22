"""Synthetic multi-segment ensembles for black-box probing of Stage-A domain selection.

Model 1 of a base PDB is cloned N times; each listed residue segment is rotated
rigidly about its own centroid by a random axis and an angle drawn from
[0.4*deg, 1.2*deg] per conformer (deg = 0 keeps the segment fixed). Internal
torsions of a segment therefore keep S = 1; only torsions spanning a segment
junction vary, so core membership is essentially every residue and the variance
clustering sees exactly the designed rigid bodies.

Optional per-segment jitter `j` adds Gaussian noise (sigma j A) to every atom of
the segment in every conformer (makes a body internally non-rigid).

Spec: comma-separated `lo-hi:deg[:j]`; residues outside every segment stay fixed.
Segments may overlap: they are applied in spec order, each about its current centroid,
so `26-76:20,51-76:5` moves 26..76 as one body and 51..76 further relative to it.
A 4th field `p<res>` pivots the rotation on the N atom of residue <res> instead of the
centroid (`39-76:6:0.05:p39` hinges at N39, so the C38-N39 peptide bond stays intact;
centroid rotation breaks the chain at every segment junction).

usage:
  segments.py gen  <out.pdb> <spec> [nconf] [seed]
  segments.py run  <name=spec> ...     # generate into cywork/seg/, run the oracle
                                         # (default, frozen, core count), print results
  segments.py runpdb <dir/x.pdb> ...   # same oracle runs on existing PDBs (one dir)
"""
import math, os, random, re, shlex, subprocess, sys
import numpy as np

BASE = "testdata/final.pdb"
ORACLE = os.environ.get("RESUPER_CYRANGE_HOST", "")  # ssh host with a licensed CYRANGE binary
REMOTE = "/tmp/resuper_seg"
OUTDIR = "cywork/seg"
RUNS = {
    "default": "",
    "frozen": "-a 1000.0 -d 50.0 -gw 0",
    "core": "-m 99 -a 1000.0 -d 50.0 -gw 0",
}


def read_model1(path):
    atoms = []
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            atoms.append((int(line[22:26]), line,
                          np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return atoms


def rot_matrix(axis, ang):
    x, y, z = axis / np.linalg.norm(axis)
    c, s = math.cos(ang), math.sin(ang)
    C = 1 - c
    return np.array([[x*x*C+c, x*y*C-z*s, x*z*C+y*s],
                     [y*x*C+z*s, y*y*C+c, y*z*C-x*s],
                     [z*x*C-y*s, z*y*C+x*s, z*z*C+c]])


def parse_spec(spec):
    """[(lo, hi, deg, jitter)] and the hinge residue per segment (None: centroid)."""
    segs, hinge = [], []
    for part in spec.split(","):
        f = part.split(":")
        lo, hi = map(int, f[0].split("-"))
        segs.append((lo, hi, float(f[1]), float(f[2]) if len(f) > 2 and f[2] else 0.0))
        hinge.append(int(f[3][1:]) if len(f) > 3 else None)
    return segs, hinge


def generate(out, spec, nconf=10, seed=7):
    rnd = random.Random(seed)
    nrnd = np.random.default_rng(seed)
    atoms = read_model1(BASE)
    segs, hinge = parse_spec(spec)
    base = np.array([c for _, _, c in atoms])
    rn = np.array([r for r, _, _ in atoms])
    names = [l[12:16].strip() for _, l, _ in atoms]
    hinge_atom = [None if h is None else
                  next(i for i in range(len(atoms)) if rn[i] == h and names[i] == "N")
                  for h in hinge]
    members = [np.where((rn >= lo) & (rn <= hi))[0] for lo, hi, _, _ in segs]
    sigma = np.zeros(len(atoms))          # jitter of the last segment holding the atom
    for (lo, hi, _, j), idx in zip(segs, members):
        sigma[idx] = j
    with open(out, "w") as f:
        f.write(f"TITLE   synthetic segments {spec}\n")
        for m in range(nconf):
            f.write(f"MODEL     {m+1:4d}\n")
            R = []
            for lo, hi, deg, _ in segs:
                ax = np.array([rnd.gauss(0, 1) for _ in range(3)])
                R.append(rot_matrix(ax, math.radians(rnd.uniform(0.4 * deg, 1.2 * deg))))
            # Segments rotate in spec order about their current centroid, so a later
            # segment nested inside an earlier one moves relative to it.
            X = base.copy()
            for Rk, idx, ha in zip(R, members, hinge_atom):
                piv = X[idx].mean(axis=0) if ha is None else X[ha].copy()
                X[idx] = (X[idx] - piv) @ Rk.T + piv
            for i, (_, l, _) in enumerate(atoms):
                p = X[i]
                if sigma[i] > 0:
                    p = p + nrnd.normal(0, sigma[i], 3)
                f.write(l[:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:])
            f.write("ENDMDL\n")
        f.write("END\n")


def summarize(text):
    out = []
    for line in text.splitlines():
        m = re.match(r"\s+Optimal range (.*)", line)
        if m:
            out.append(m.group(1))
        m = re.search(r"number of core atoms \(\s*(\d+)\)", line)
        if m:
            out.append(f"C={m.group(1)}")
        if "Only domains" in line or "egmentation" in line:
            out.append(line.strip())
    return " | ".join(out) if out else text.strip()[-200:]


def run(pairs):
    os.makedirs(OUTDIR, exist_ok=True)
    names = []
    for pair in pairs:
        name, spec = pair.split("=", 1)
        generate(os.path.join(OUTDIR, name + ".pdb"), spec)
        names.append((name, spec))
    oracle(names, OUTDIR)


def oracle(names, outdir):
    """Run the three oracle modes on <outdir>/<name>.pdb for every (name, label)."""
    OUTDIR = outdir
    subprocess.run(["ssh", ORACLE, f"mkdir -p {REMOTE} && rm -f {REMOTE}/*"], check=True)
    subprocess.run(["scp", "-q"] + [os.path.join(OUTDIR, n + ".pdb") for n, _ in names]
                   + [f"{ORACLE}:{REMOTE}/"], check=True)
    body = f"cd {REMOTE} && for p in {' '.join(n for n, _ in names)}; do "
    for tag, flags in RUNS.items():
        body += f"/usr/local/bin/cyrange {flags} $p.pdb > $p.{tag} 2>&1; "
    body += "done"
    subprocess.run(["ssh", ORACLE, body], check=True)
    tar = subprocess.run(["ssh", ORACLE, f"cd {REMOTE} && tar cf - --exclude='*.pdb' ."],
                         capture_output=True, check=True).stdout
    subprocess.run(["tar", "xf", "-", "-C", OUTDIR], input=tar, check=True)
    for name, spec in names:
        print(f"== {name}  {spec}")
        for tag in RUNS:
            p = os.path.join(OUTDIR, f"{name}.{tag}")
            print(f"   {tag:7s} {summarize(open(p).read()) if os.path.exists(p) else 'MISSING'}")


if __name__ == "__main__":
    if sys.argv[1] == "gen":
        a = sys.argv[2:]
        generate(a[0], a[1], int(a[2]) if len(a) > 2 else 10, int(a[3]) if len(a) > 3 else 7)
    elif sys.argv[1] == "run":
        run(sys.argv[2:])
    elif sys.argv[1] == "runpdb":
        paths = sys.argv[2:]
        oracle([(os.path.basename(p)[:-4], p) for p in paths], os.path.dirname(paths[0]))
