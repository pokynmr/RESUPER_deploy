"""Fetch independent domain annotations (CATH, SCOP, ECOD) for PDB entries.

These classifications are made independently of RESUPER and CYRANGE, so they can
serve as independent evidence that a protein has more than one structural
domain. Caveats:
  * classification domains are folding/evolutionary units, not necessarily
    rigid bodies that move independently (e.g. 1DX7: paper describes two
    helical domains, all three classifications annotate one);
  * positions are in the PDB's sequential numbering (label_seq_id, starting at 1);
    this script also records the mapping to the author residue numbers used in
    the coordinate files and in RESUPER's reports.

Source: RCSB Data API, https://data.rcsb.org/rest/v1/core/polymer_entity_instance/<ID>/<chain>
Raw responses are saved, so the table can be re-derived without the network.

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/domain_annotations.py OUTDIR ID[:CHAIN] ...
Output: OUTDIR/raw/<ID>_<CHAIN>.json, OUTDIR/domain_annotations.csv
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True

import bench_config as cfg                       # noqa: E402
from fetch_structures import http_get            # noqa: E402

URL = "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{id}/{chain}"
TYPES = ("CATH", "SCOP", "ECOD")


def author_map(d):
    """label_seq_id -> author residue number, if the API provides it."""
    m = {}
    ids = d.get("rcsb_polymer_entity_instance_container_identifiers", {})
    seq = ids.get("auth_to_entity_poly_seq_mapping")
    if seq:
        for i, a in enumerate(seq, start=1):
            m[i] = a
    return m


def main():
    log = cfg.setup_logging("domain_annotations")
    outdir = sys.argv[1]
    os.makedirs(os.path.join(outdir, "raw"), exist_ok=True)
    rows = []
    for arg in sys.argv[2:]:
        pid, chain = (arg.split(":") + ["A"])[:2]
        raw = os.path.join(outdir, "raw", f"{pid}_{chain}.json")
        try:
            if not os.path.exists(raw):
                data = http_get(URL.format(id=pid, chain=chain), log)
                json.loads(data)
                open(raw, "wb").write(data)
            d = json.load(open(raw))
        except Exception as e:
            rows.append({"pdb_id": pid, "chain": chain, "status": f"failed: {type(e).__name__}: {e}"})
            log.error("%s %s: %s", pid, chain, e)
            continue
        amap = author_map(d)
        row = {"pdb_id": pid, "chain": chain, "status": "ok",
               "author_numbering_available": bool(amap)}
        for t in TYPES:
            # One feature per segment; segments of the same domain share a
            # *_DOMAIN_ID property (e.g. CATH 2d21A01), so group by it.
            doms = {}
            for n, f in enumerate(d.get("rcsb_polymer_instance_feature", [])):
                if f["type"] != t:
                    continue
                dom_id = next((v["values"][0] for v in f.get("additional_properties") or []
                               if v["name"].endswith("DOMAIN_ID") and v.get("values")), f"feature{n}")
                entry = doms.setdefault(dom_id, {"name": f.get("name") or "", "segs": []})
                for p in f.get("feature_positions", []):
                    a, b = p["beg_seq_id"], p["end_seq_id"]
                    entry["segs"].append(f"{a}-{b} (author {amap.get(a, '?')}-{amap.get(b, '?')})"
                                         if amap else f"{a}-{b}")
            row[f"{t}_n_domains"] = len(doms)
            row[f"{t}_domains"] = " ; ".join(f"{k} {v['name']} [{' + '.join(v['segs'])}]"
                                             for k, v in doms.items())
        rows.append(row)
        log.info("%s %s: CATH %s, SCOP %s, ECOD %s", pid, chain, row["CATH_n_domains"],
                 row["SCOP_n_domains"], row["ECOD_n_domains"])
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(os.path.join(outdir, "domain_annotations.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
