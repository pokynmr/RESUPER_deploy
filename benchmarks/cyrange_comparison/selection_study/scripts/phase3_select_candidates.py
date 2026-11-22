"""Selection study, Phase 3: choose candidate NMR structures by independent properties.

Implements decision S6 (docs/decisions_log.md), fixed before this script was run.
Nothing here downloads coordinates or runs RESUPER: candidates are chosen only
from PDB metadata, so the choice cannot depend on how either method performs.

Steps
  1. RCSB Search API: all released SOLUTION NMR entries with >= 10 models, one
     protein entity, one polymer instance, no nucleic acid.
  2. Remove every PDB ID used in RESUPER's development (repository test data and
     benchmark lists) and the 37 CYRANGE-paper proteins analysed earlier.
  3. RCSB GraphQL Data API (batches of 100): title, organism, length, model
     count, primary citation, CATH/SCOP/ECOD domains, helices and strands.
     Raw responses are cached (results/phase3/raw/) so the selection can be
     re-derived offline.
  4. Assign each entry to at most one group (B, A, D, C, E; rules in S6).
  5. Draw 12 per group at random (numpy seed 20261111).
  6. For selected entries with a PubMed ID, fetch the abstract from Europe PMC and
     record sentences mentioning flexibility-related terms (leads only).
  7. Write the pool, the selection, and the selection's SHA-256 (frozen list).

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/phase3_select_candidates.py
"""
import csv
import glob
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True

import numpy as np                    # noqa: E402
import bench_config as cfg            # noqa: E402

OUT = os.path.join(cfg.BENCH_DIR, "selection_study", "results", "phase3")
RAW = os.path.join(OUT, "raw")
SEED, PER_GROUP = 20261111, 12
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
TERMS = re.compile(r"flexib|hinge|inter-?domain|linker|dynamic|relaxation|motion|mobil|"
                   r"orientation|tumbl|independen|disorder", re.I)

GQL = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
  struct{title}
  rcsb_entry_info{deposited_model_count deposited_polymer_monomer_count}
  rcsb_accession_info{initial_release_date}
  rcsb_primary_citation{title journal_abbrev year pdbx_database_id_PubMed pdbx_database_id_DOI}
  polymer_entities{ rcsb_entity_source_organism{scientific_name}
    rcsb_polymer_entity{pdbx_description} entity_poly{rcsb_sample_sequence_length}
    polymer_entity_instances{ rcsb_polymer_entity_instance_container_identifiers{auth_asym_id}
      rcsb_polymer_instance_feature{ type name feature_positions{beg_seq_id end_seq_id}
        additional_properties{name values} } } } } }"""


def post(url, payload, log, tries=4):
    data = json.dumps(payload).encode()
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={
                "Content-Type": "application/json", "User-Agent": cfg.USER_AGENT})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except Exception as e:                       # network: retry with backoff
            if k == tries - 1:
                raise
            log.warning("%s: %s; retry in %d s", url, e, 2 ** (k + 1))
            time.sleep(2 ** (k + 1))


def get(url, log, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def search_pool(log):
    def t(attr, op, val):
        return {"type": "terminal", "service": "text",
                "parameters": {"attribute": attr, "operator": op, "value": val}}
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        t("exptl.method", "exact_match", "SOLUTION NMR"),
        t("rcsb_entry_info.deposited_model_count", "greater_or_equal", 10),
        t("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
        t("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
        t("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0)]},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    res = post(SEARCH, q, log)
    json.dump({"query": q, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
               "total_count": res["total_count"]}, open(os.path.join(RAW, "search_query.json"), "w"), indent=1)
    return sorted(r["identifier"] for r in res["result_set"])


def excluded_ids():
    """PDB IDs from RESUPER's development material + the earlier 37-protein benchmark."""
    ids = set()
    pat = re.compile(r"^([0-9][A-Za-z0-9]{3})(?:[_.]|$)")
    for f in glob.glob(os.path.join(cfg.REPO_TESTDATA, "**", "*"), recursive=True):
        m = pat.match(os.path.basename(f))
        if m:
            ids.add(m.group(1).upper())
    for lst in (os.path.join(cfg.REPO_TESTDATA, "pacsy_ids.txt"),
                os.path.join(cfg.REPO_DIR, "scripts", "bench", "proteins.txt")):
        if os.path.exists(lst):
            for line in open(lst):
                m = pat.match(line.strip())
                if m:
                    ids.add(m.group(1).upper())
    for row in csv.DictReader(open(cfg.MANIFEST)):
        ids.add(row["pdb_id"].upper())
    return ids


def fetch_metadata(ids, log):
    out = {}
    for i in range(0, len(ids), 100):
        batch = ids[i:i + 100]
        cache = os.path.join(RAW, f"graphql_{i // 100:04d}.json.gz")
        if os.path.exists(cache):
            data = json.load(gzip.open(cache, "rt"))
        else:
            data = post(GRAPHQL, {"query": GQL, "variables": {"ids": batch}}, log)
            json.dump(data, gzip.open(cache, "wt"))
        for e in (data.get("data") or {}).get("entries") or []:
            if e:
                out[e["rcsb_id"]] = e
        if (i // 100) % 20 == 0:
            log.info("metadata: %d / %d", min(i + 100, len(ids)), len(ids))
    return out


def summarise(e):
    pe = e["polymer_entities"][0]
    inst = pe["polymer_entity_instances"][0]
    feats = inst["rcsb_polymer_instance_feature"] or []
    L = pe["entity_poly"]["rcsb_sample_sequence_length"]
    doms = {}
    for t in ("CATH", "SCOP", "ECOD"):
        d = {}
        for n, f in enumerate(feats):
            if f["type"] != t:
                continue
            key = next((v["values"][0] for v in f.get("additional_properties") or []
                        if v["name"].endswith("DOMAIN_ID") and v.get("values")), f"f{n}")
            d.setdefault(key, set()).update(
                r for p in f["feature_positions"] for r in range(p["beg_seq_id"], p["end_seq_id"] + 1))
        doms[t] = list(d.values())
    helices = sorted((p["beg_seq_id"], p["end_seq_id"]) for f in feats if f["type"] == "HELIX_P"
                     for p in f["feature_positions"])
    strands = sorted((p["beg_seq_id"], p["end_seq_id"]) for f in feats if f["type"] == "SHEET"
                     for p in f["feature_positions"])
    cit = e.get("rcsb_primary_citation") or {}
    org = pe.get("rcsb_entity_source_organism") or []
    return {"pdb_id": e["rcsb_id"], "title": (e.get("struct") or {}).get("title", ""),
            "protein": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", ""),
            "organism": "; ".join(sorted({o["scientific_name"] for o in org if o.get("scientific_name")})),
            "length": L, "n_models": e["rcsb_entry_info"]["deposited_model_count"],
            "chain": (inst.get("rcsb_polymer_entity_instance_container_identifiers") or {}).get("auth_asym_id", ""),
            "released": (e.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
            "citation": f"{cit.get('title') or ''} {cit.get('journal_abbrev') or ''} {cit.get('year') or ''}".strip(),
            "pmid": cit.get("pdbx_database_id_PubMed") or "", "doi": cit.get("pdbx_database_id_DOI") or "",
            "_doms": doms, "_helices": helices, "_strands": strands}


def assign(s):
    L = s["length"]
    counts = {t: len(v) for t, v in s["_doms"].items() if v}
    s["n_dom_CATH"], s["n_dom_SCOP"], s["n_dom_ECOD"] = (len(s["_doms"][t]) for t in ("CATH", "SCOP", "ECOD"))
    union = set().union(*[r for v in s["_doms"].values() for r in v]) if counts else set()
    s["residues_outside_domains"] = L - len(union & set(range(1, L + 1))) if counts else ""
    s["helices"] = " ".join(f"{a}-{b}" for a, b in s["_helices"])
    s["n_strands"] = len(s["_strands"])
    vals = list(counts.values())
    multi = len(vals) >= 2 and min(vals) >= 2
    med = float(np.median(vals)) if vals else 0
    if multi and med >= 3 and 120 <= L <= 500:
        return "B", f"{len(vals)} sources agree on >= 2 domains, median {med:g}"
    if multi and med == 2 and 80 <= L <= 500:
        return "A", f"{len(vals)} sources agree on 2 domains"
    if vals and max(vals) == 1 and 60 <= L <= 400:
        for a, b in s["_helices"]:
            n = b - a + 1
            if 8 <= n <= 25 and all(r not in union for r in range(a - 3, b + 4)):
                return "D", f"single domain; helix {a}-{b} ({n} res) outside it"
    if len(vals) >= 2 and max(vals) == 1 and 60 <= L <= 300:
        out = s["residues_outside_domains"]
        if out >= max(20, 0.2 * L):
            return "C", f"single domain; {out} of {L} residues outside every annotated domain"
    if 25 <= L <= 50 and s["_helices"]:
        longest = max(b - a + 1 for a, b in s["_helices"])
        if longest >= 0.7 * L:
            return "E", f"peptide; longest helix {longest} of {L} residues"
    return "", ""


def literature(pmid, log):
    if not pmid:
        return "", ""
    try:
        q = urllib.parse.urlencode({"query": f"EXT_ID:{pmid} AND SRC:MED", "format": "json",
                                    "resultType": "core"})
        r = get(f"{EPMC}?{q}", log)
        res = (r.get("resultList") or {}).get("result") or []
        abstract = re.sub(r"<[^>]+>", "", (res[0].get("abstractText") or "")) if res else ""
    except Exception as e:
        return "", f"abstract unavailable: {type(e).__name__}"
    sents = [s.strip() for s in re.split(r"(?<=[.!?])\s+", abstract) if TERMS.search(s)]
    return ("yes" if abstract else "no abstract"), " | ".join(sents[:3])


def main():
    if os.path.exists(os.path.join(OUT, "FROZEN_selection.txt")) and "--redraw" not in sys.argv:
        print("Frozen selection exists (FROZEN_selection.txt); not redrawing. The PDB changes over time, so a redraw "
              "would select different entries. Pass --redraw to draw a new selection.")
        return 0
    log = cfg.setup_logging("selection_phase3")
    os.makedirs(RAW, exist_ok=True)
    pool = search_pool(log)
    excl = excluded_ids()
    ids = [i for i in pool if i not in excl]
    log.info("search: %d entries; %d excluded as RESUPER development or earlier benchmark; %d remain",
             len(pool), len(pool) - len(ids), len(ids))
    meta = fetch_metadata(ids, log)
    rows = []
    for pid in ids:
        if pid not in meta:
            continue
        s = summarise(meta[pid])
        s["group"], s["group_reason"] = assign(s)
        rows.append(s)
    pub = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    fields = list(pub[0].keys())
    with open(os.path.join(OUT, "candidate_pool.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(pub)
    sizes = {g: sum(r["group"] == g for r in rows) for g in "ABCDE"}
    log.info("eligible pool sizes: %s", sizes)

    rng = np.random.default_rng(SEED)
    selected = []
    for g in "BADCE":
        members = sorted((r for r in pub if r["group"] == g), key=lambda r: r["pdb_id"])
        k = min(PER_GROUP, len(members))
        pick = rng.choice(len(members), size=k, replace=False) if k else []
        for i in sorted(pick):
            selected.append(dict(members[i]))
    for r in selected:
        r["abstract_found"], r["abstract_flexibility_sentences"] = literature(r["pmid"], log)
    order = ["group", "pdb_id", "protein", "organism", "length", "n_models", "chain", "group_reason",
             "n_dom_CATH", "n_dom_SCOP", "n_dom_ECOD", "residues_outside_domains", "helices",
             "n_strands", "title", "released", "citation", "pmid", "doi", "abstract_found",
             "abstract_flexibility_sentences"]
    path = os.path.join(OUT, "candidates_selected.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=order, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(selected, key=lambda r: ("BADCE".index(r["group"]), r["pdb_id"])))
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(os.path.join(OUT, "FROZEN_selection.txt"), "w") as fh:
        fh.write(f"candidates_selected.csv sha256 {digest}\nfrozen {datetime.now(timezone.utc).isoformat()}\n"
                 f"pool sizes {sizes}; seed {SEED}; {PER_GROUP} per group\n")
    log.info("selected %d candidates; frozen sha256 %s", len(selected), digest[:16])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
