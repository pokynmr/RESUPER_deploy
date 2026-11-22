"""Select a diverse set of solution-NMR multi-model proteins directly from the PDB.

Queries the public RCSB PDB search and data APIs (no local database, SSH or
credentials): NMR entries with at least two deposited models and a single
polymer chain instance. The candidate pool is then deterministically sampled
across (residue-count, model-count) buckets, preferring larger and multi-model
entries so that edge cases are exercised.

The committed testdata/pacsy_ids.txt is a frozen historical benchmark selected
from PACSY metadata in 2026. RCSB cannot reproduce that exact set (four of its
entries are absent from the current single-chain pool), so this script is for
new selections and does not rewrite the frozen list.

usage: /usr/bin/python3 scripts/select_pdb.py [--count K] [--out PATH]
  default: print the selected ids to stdout
  --out PATH: also write them, one per line
"""
import json
import os
import sys
import urllib.request

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
PAGE = 10000
CHUNK = 500

QUERY = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.experimental_method",
                "operator": "exact_match", "value": "NMR"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.deposited_model_count",
                "operator": "greater_or_equal", "value": 2}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute":
                    "rcsb_entry_info.deposited_polymer_entity_instance_count",
                "operator": "equals", "value": 1}},
        ],
    },
    "return_type": "entry",
}


def _post(url, payload):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.load(r)


def fetch_ids():
    ids, start = [], 0
    while True:
        body = dict(QUERY)
        body["request_options"] = {"paginate": {"start": start, "rows": PAGE}}
        d = _post(SEARCH, body)
        ids += [x["identifier"] for x in d["result_set"]]
        if len(ids) >= d["total_count"] or not d["result_set"]:
            break
        start += len(d["result_set"])
    return ids


def fetch_metadata(ids):
    out = {}
    for i in range(0, len(ids), CHUNK):
        chunk = ids[i:i + CHUNK]
        q = ("{ entries(entry_ids: ["
             + ",".join('"%s"' % x for x in chunk)
             + "]) { rcsb_id rcsb_entry_info { deposited_model_count "
               "deposited_polymer_monomer_count "
               "deposited_polymer_entity_instance_count "
               "experimental_method } } }")
        d = _post(GRAPHQL, {"query": q})
        if "errors" in d:
            raise RuntimeError(d["errors"])
        for e in d["data"]["entries"]:
            info = e["rcsb_entry_info"]
            out[e["rcsb_id"]] = (info["deposited_model_count"],
                                 info["deposited_polymer_monomer_count"],
                                 info["deposited_polymer_entity_instance_count"])
    return out


def select(rep, k):
    ids = sorted(rep)

    def sizeb(sc):
        return 0 if sc < 50 else (1 if sc < 100 else (2 if sc < 160 else 3))

    def modelb(mc):
        return 0 if mc <= 10 else (1 if mc <= 25 else 2)

    buckets = {}
    for pid in ids:
        mc, sc, nc = rep[pid]
        buckets.setdefault((sizeb(sc), modelb(mc), 1 if nc >= 2 else 0),
                           []).append(pid)

    chosen = []
    taglist = sorted(buckets.keys())
    taglist.sort(key=lambda t: (-1 if t[2] else 0, -t[0], t[1]))
    bytag = {t: sorted(buckets[t]) for t in taglist}
    idx = {t: 0 for t in taglist}
    active = list(t for t in taglist if bytag[t])
    while len(chosen) < k and active:
        nxt = []
        for t in list(active):
            if len(chosen) >= k:
                break
            if idx[t] < len(bytag[t]):
                chosen.append(bytag[t][idx[t]])
                idx[t] += 1
                nxt.append(t)
        active = [t for t in nxt if idx[t] < len(bytag[t])]
    return chosen


def main():
    k = 100
    out = None
    args = sys.argv[1:]
    if "--count" in args:
        k = int(args[args.index("--count") + 1])
    if "--out" in args:
        out = args[args.index("--out") + 1]

    rep = fetch_metadata(fetch_ids())
    chosen = select(rep, k)
    print(f"selected {len(chosen)} ids from {len(rep)} candidates",
          file=sys.stderr)
    for pid in chosen:
        print(pid)
    if out:
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        with open(out, "w") as f:
            f.write("\n".join(chosen) + "\n")
        print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
