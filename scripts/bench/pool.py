"""Protein pool for benchmark B: solution-NMR monomers sampled from the PDB.

Query (RCSB search API): solution NMR, one protein entity, one polymer instance, no
nucleic acid, 80-250 deposited residues, at least 10 models. The hits are sorted by ID,
shuffled with a fixed seed, and taken in order, skipping IDs already in testdata/, until
POOL_SIZE entries pass the local checks: one chain, standard residues only, no residue-
number gaps, every model with the same atoms, at least 10 models. The chosen IDs are
written to scripts/bench/proteins.txt (committed); the files are cached in cywork/bench/pdb/.

usage: /usr/bin/python3 scripts/bench/pool.py [--seed N --list FILE]   # select, download
       /usr/bin/python3 scripts/bench/pool.py --fetch [--list FILE]      # re-download
"""
import json, os, random, sys, urllib.request

POOL_SIZE = 40
SEED = 20261021
CACHE = "cywork/bench/pdb"
LIST = "scripts/bench/proteins.txt"
# Benchmark C: --seed 20261022 --list scripts/bench/proteins_c.txt, excluding the IDs of
# every other list in scripts/bench/ as well as testdata/.
STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
            "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def query():
    t = lambda attr, op, val: {"type": "terminal", "service": "text",
                               "parameters": {"attribute": attr, "operator": op, "value": val}}
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        t("exptl.method", "exact_match", "SOLUTION NMR"),
        t("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
        t("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
        t("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
        t("rcsb_entry_info.deposited_polymer_monomer_count", "range",
          {"from": 80, "to": 250, "include_lower": True, "include_upper": True}),
        t("rcsb_entry_info.deposited_model_count", "greater_or_equal", 10)]},
         "return_type": "entry", "request_options": {"return_all_hits": True}}
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
                                 data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return sorted(h["identifier"] for h in json.load(r)["result_set"])


def download(pid):
    path = f"{CACHE}/{pid}.pdb"
    if not os.path.exists(path):
        req = urllib.request.Request(f"https://files.rcsb.org/download/{pid}.pdb",
                                     headers={"User-Agent": "resuper"})
        with urllib.request.urlopen(req, timeout=120) as r:
            open(path, "wb").write(r.read())
    return path


def usable(path):
    models, cur = [], None
    for line in open(path):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = None
        elif cur is not None and line.startswith(("ATOM", "HETATM")):
            cur.append(line)
    if len(models) < 10:
        return False
    first = models[0]
    if any(l.startswith("HETATM") and l[17:20] != "HOH" for l in first):
        return False
    atoms = [l for l in first if l.startswith("ATOM")]
    if len({l[21] for l in atoms}) != 1 or any(l[16] not in " " for l in atoms):
        return False
    if any(l[17:20] not in STANDARD for l in atoms) or any(l[26] != " " for l in atoms):
        return False
    nums = [int(l[22:26]) for l in atoms if l[12:16] == " CA "]
    if any(b != a + 1 for a, b in zip(nums, nums[1:])):
        return False
    key = [l[12:27] for l in atoms]
    return all([l[12:27] for l in m if l.startswith("ATOM")] == key for m in models)


def main():
    global SEED, LIST
    if "--seed" in sys.argv:
        SEED = int(sys.argv[sys.argv.index("--seed") + 1])
    if "--list" in sys.argv:
        LIST = sys.argv[sys.argv.index("--list") + 1]
    os.makedirs(CACHE, exist_ok=True)
    if "--fetch" in sys.argv:
        for pid in open(LIST).read().split():
            download(pid)
        return
    have = {os.path.basename(f)[:-4].upper() for f in os.listdir("testdata") if f.endswith(".pdb")}
    for f in os.listdir("scripts/bench"):
        if f.startswith("proteins") and f.endswith(".txt") and f"scripts/bench/{f}" != LIST:
            have |= {x.upper() for x in open(f"scripts/bench/{f}").read().split()}
    ids = query()
    random.Random(SEED).shuffle(ids)
    chosen = []
    for pid in ids:
        if pid.upper() in have:
            continue
        try:
            if usable(download(pid)):
                chosen.append(pid)
        except Exception as e:
            print(f"skip {pid}: {e}")
        if len(chosen) == POOL_SIZE:
            break
    open(LIST, "w").write("\n".join(chosen) + "\n")
    print(f"{len(ids)} hits; chose {len(chosen)}: {' '.join(chosen)}")


if __name__ == "__main__":
    main()
