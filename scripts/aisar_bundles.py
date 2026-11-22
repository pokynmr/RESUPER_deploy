"""Build multi-model PDB bundles from the AISAR two-state ensembles (Huang et al.,
Nat Commun 2026; https://github.com/MontelioneLab/AISAR) for oracle Stage-A probing.

Inputs are fetched into cywork/aisar/raw/ (repo-relative paths kept, pinned to
commit REV). Each bundle concatenates the ATOM records of every model of every listed
source file, renumbered as MODEL 1..N, so a two-state mixture looks like one NMR
ensemble with a large lid motion.

The bundles are held-out regression cases (CASES in scripts/eval_cases.py and
tests/test_reference.py): their oracle reports are committed as testdata/<name>.out,
the PDBs are not (testdata/.gitignore).

usage: /usr/bin/python3 scripts/aisar_bundles.py [--fetch] [--testdata]
  --fetch     download the source files into cywork/aisar/raw/
  --testdata  also copy the bundles into testdata/
"""
import glob, os, shutil, sys, urllib.request

REV = "13d0ba89b84c8e87168c451defab390525d77ae1"
URL = f"https://raw.githubusercontent.com/MontelioneLab/AISAR/{REV}/"

RAW = "cywork/aisar/raw"
OUT = "cywork/aisar"
G = "GLuc/figs_tables/selectedModels-20"
C = "CDK2AP1/figs_tables/selectedModels"
DR = "CDK2AP1/doubleRecall/TableS7_refinedNOESYPeakList"
BUNDLES = {
    "gluc_s1": ["GLuc/doubleRecall/class1.pdb"],
    "gluc_s2": ["GLuc/doubleRecall/class2.pdb"],
    "gluc_s12": ["GLuc/doubleRecall/class1.pdb", "GLuc/doubleRecall/class2.pdb"],
    "gluc_nmr": ["GLuc/doubleRecall/nmr.pdb"],
    "gluc_a5": [f"{G}/a5/*.pdb"],
    "gluc_b3": [f"{G}/b3/*.pdb"],
    "gluc_a5b3": [f"{G}/a5/*.pdb", f"{G}/b3/*.pdb"],
    "gluc_c2": [f"{G}/c2/*.pdb"],
    "gluc_d4": [f"{G}/d4/*.pdb"],
    "cdk_s1": [f"{DR}/state1.pdb"],
    "cdk_s2": [f"{DR}/state2.pdb"],
    "cdk_s12": [f"{DR}/state1.pdb", f"{DR}/state2.pdb"],
    "cdk_nmr": [f"{DR}/nmr.pdb"],
    "cdk_a1b2": [f"{C}/a1/*.pdb", f"{C}/b2/*.pdb"],
    "cdk_c4d3": [f"{C}/c4/*.pdb", f"{C}/d3/*.pdb"],
}


def sources():
    """Repo paths of every source file the bundles use."""
    out = []
    for pats in BUNDLES.values():
        for pat in pats:
            if "*" not in pat:
                out.append(pat)
            elif pat.startswith(G):
                out += [pat.replace("*", f"{pat.split('/')[-2]}_{i}") for i in range(10, 21)]
            else:
                out += [pat.replace("*", f"{pat.split('/')[-2]}_{i}") for i in range(1, 6)]
    return sorted(set(out))


def fetch():
    for rel in sources():
        dst = os.path.join(RAW, rel)
        if os.path.exists(dst):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with urllib.request.urlopen(URL + rel) as r, open(dst, "wb") as fh:
            fh.write(r.read())


def models(path):
    """ATOM blocks of each model in a PDB file (a file without MODEL is one model)."""
    out, cur = [], []
    for line in open(path):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ATOM"):
            cur.append(line if line.endswith("\n") else line + "\n")
        elif line.startswith("ENDMDL"):
            out.append(cur)
            cur = []
    if cur:
        out.append(cur)
    return out


def main():
    if "--fetch" in sys.argv:
        fetch()
    for name, pats in BUNDLES.items():
        blocks = []
        for pat in pats:
            for f in sorted(glob.glob(os.path.join(RAW, pat))):
                blocks.extend(models(f))
        sizes = {len(b) for b in blocks}
        with open(os.path.join(OUT, name + ".pdb"), "w") as fh:
            for k, b in enumerate(blocks, 1):
                fh.write(f"MODEL     {k:4d}\n")
                fh.writelines(b)
                fh.write("ENDMDL\n")
            fh.write("END\n")
        print(f"{name:10s} models={len(blocks):3d} atoms/model={sorted(sizes)}")
        if "--testdata" in sys.argv:
            shutil.copy(os.path.join(OUT, name + ".pdb"), os.path.join("testdata", name + ".pdb"))


if __name__ == "__main__":
    main()
