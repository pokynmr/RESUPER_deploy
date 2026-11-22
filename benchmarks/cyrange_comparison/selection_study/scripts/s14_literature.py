"""Literature for the top-ranked expanded-scan cases (S14).
For the first N entries of results/expanded/ranked.csv: Europe PMC abstract by PMID
(from the Phase 3 pool) or by PDB ID, plus open-access full-text paragraphs that match
flexibility/domain/helix terms. Saved unmodified as text for reading.
usage (from benchmarks/cyrange_comparison): python selection_study/scripts/s14_literature.py [N]
Output: results/expanded/literature/<ID>.txt
"""
import csv, json, os, re, sys, urllib.parse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
from phase6_literature import EPMC, TERMS, fetch, paragraphs  # noqa: E402
import bench_config as cfg  # noqa: E402
S = os.path.join(cfg.BENCH_DIR, "selection_study", "results")
OUT = os.path.join(S, "expanded", "literature")
EXTRA = re.compile(r"helix|helical|kink|N-termin|C-termin|unstructured|two domains|three domains|"
                   r"lobe|subdomain|domain", re.I)
# argument: N (top N of ranked.csv) or a comma-separated list of PDB IDs
ARG = sys.argv[1] if len(sys.argv) > 1 else "30"
IDS = None if ARG.isdigit() else set(ARG.upper().split(","))
N = int(ARG) if IDS is None else None
pool = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(S, "phase3", "candidate_pool.csv")))}
os.makedirs(OUT, exist_ok=True)
for r in [x for x in csv.DictReader(open(os.path.join(S, "expanded", "ranked.csv"))) if IDS is None or x["pdb_id"] in IDS][:N]:
    pid = r["pdb_id"]; p = pool.get(pid, {})
    q = f"EXT_ID:{p['pmid']}" if p.get("pmid") else f'"{pid}"'
    lines = [f"{pid}  pool citation: {p.get('citation', '')[:200]}  query: {q}"]
    try:
        res = json.loads(fetch(f"{EPMC}/search?query={urllib.parse.quote(q)}&resultType=core&format=json&pageSize=1"))
        hits = res.get("resultList", {}).get("result", [])
    except Exception as e:  # noqa: BLE001
        hits = []; lines.append(f"search failed: {e}")
    for h in hits[:1]:
        lines += [f"TITLE: {h.get('title')}", f"YEAR {h.get('pubYear')}  PMID {h.get('pmid')}  PMCID {h.get('pmcid')}",
                  "ABSTRACT: " + re.sub(r"<[^>]+>", "", h.get("abstractText", "") or "(none)")]
        if h.get("pmcid") and h.get("isOpenAccess") == "Y":
            try:
                paras = paragraphs(fetch(f"{EPMC}/{h['pmcid']}/fullTextXML"))
                sel = [(t, x) for t, x in paras if TERMS.search(x) or EXTRA.search(x)]
                lines += [f"FULL TEXT: {len(paras)} paragraphs, {len(sel)} matching"] + [f"[{t}] {x}" for t, x in sel]
            except Exception as e:  # noqa: BLE001
                lines.append(f"full text failed: {e}")
    if not hits:
        lines.append("(no paper found)")
    open(os.path.join(OUT, f"{pid}.txt"), "w").write("\n".join(lines) + "\n")
    print(pid, "paper" if hits else "none", "full text" if any(l.startswith("FULL TEXT") for l in lines) else "")
