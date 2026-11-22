"""Build a synthetic TWO-DOMAIN ensemble from final.pdb: split the chain into two
rigid halves (domain A = res 1..bnd, domain B = bnd+1..76) and, per conformer,
rotate each half independently around its own pivot. The junction becomes mobile,
so CYRANGE Stage A should identify two separate domains.

usage: /usr/bin/python3 scripts/make_2domain.py [nconf] [bnd] [deg] [out]
"""
import math, random, sys
import numpy as np

def read_atoms(path):
    lines = []
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            resnum = int(line[22:26])
            name = line[12:16].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            lines.append((resnum, line, np.array([x, y, z])))
    return lines

def rot_matrix(axis, ang):
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = math.cos(ang); s = math.sin(ang); C = 1 - c
    return np.array([
        [x*x*C+c, x*y*C-z*s, x*z*C+y*s],
        [y*x*C+z*s, y*y*C+c, y*z*C-x*s],
        [z*x*C-y*s, z*y*C+x*s, z*z*C+c]])

def main():
    nconf = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    bnd = int(sys.argv[2]) if len(sys.argv) > 2 else 38
    deg = float(sys.argv[3]) if len(sys.argv) > 3 else 35.0
    out = sys.argv[4] if len(sys.argv) > 4 else "testdata/two_domain.pdb"
    random.seed(7)
    atoms = read_atoms("testdata/final.pdb")
    A = [(rn, l, c) for (rn, l, c) in atoms if rn <= bnd]
    B = [(rn, l, c) for (rn, l, c) in atoms if rn > bnd]
    pivA = np.mean([c for _, _, c in A], axis=0)
    pivB = np.mean([c for _, _, c in B], axis=0)
    with open(out, "w") as f:
        f.write("TITLE   synthetic two-domain (rigid halves, flexible junction)\n")
        for k in range(nconf):
            f.write(f"MODEL     {k+1:4d}\n")
            # per-conformer independent random rotations for each half
            axA = np.array([random.gauss(0, 1), random.gauss(0, 1), random.gauss(0, 1)])
            axB = np.array([random.gauss(0, 1), random.gauss(0, 1), random.gauss(0, 1)])
            RA = rot_matrix(axA, math.radians(random.uniform(0.4 * deg, 1.2 * deg)))
            RB = rot_matrix(axB, math.radians(random.uniform(0.4 * deg, 1.2 * deg)))
            for rn, l, c in A:
                p = RA @ (c - pivA) + pivA
                # re-emit with new coords
                f.write(l[:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:])
            for rn, l, c in B:
                p = RB @ (c - pivB) + pivB
                f.write(l[:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:])
            f.write("ENDMDL\n")
        f.write("END\n")
    print("wrote", out, "nconf=", nconf, "bnd=", bnd, "deg=", deg)

if __name__ == "__main__":
    main()
