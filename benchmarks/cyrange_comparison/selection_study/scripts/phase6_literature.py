"""Selection study, Phase 6a: retrieve the primary papers of the Phase 6 proteins.

For each PubMed ID, Europe PMC is queried for the article record. If an
open-access full text exists (PMCID + isOpenAccess), its XML is downloaded;
otherwise the abstract is used and this is recorded. Paragraphs that mention
dynamics or domain organisation are extracted for reading. Nothing is
interpreted here: the excerpts are read and summarised in the Phase 6 write-up.

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/phase6_literature.py
Output: selection_study/results/phase6/literature/<ID>_record.json,
        <ID>_fulltext.xml (if open access), <ID>_excerpts.txt, literature_status.csv
"""
import csv
import json
import os
import re
import sys
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
OUT = os.path.join(S, "results", "phase6", "literature")
CASES = ["1OQY", "2P01", "2KUI", "2BW2", "2RPC", "2F3J", "7A2D", "2LQ4", "2M70", "1B8T"]
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
TERMS = re.compile(r"relaxation|T1\b|T2\b|heteronuclear NOE|\{1H\}|order parameter|S2\b|"
                   r"dipolar|RDC|SAXS|scattering|tumbl|correlation time|linker|inter-?domain|"
                   r"flexib|hinge|rigid|orientation|independen|mobil|dynamic|disorder", re.I)


def fetch(url, binary=False):
    req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    return data if binary else data.decode("utf-8", "replace")


def paragraphs(xml):
    """Plain-text paragraphs (with section titles) from a JATS full-text XML."""
    out = []
    for sec in re.findall(r"<sec[^>]*>(.*?)</sec>", xml, flags=re.S) or [xml]:
        title = re.search(r"<title>(.*?)</title>", sec, flags=re.S)
        t = re.sub(r"<[^>]+>", "", title.group(1)).strip() if title else ""
        for p in re.findall(r"<p[^>]*>(.*?)</p>", sec, flags=re.S):
            txt = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", p)).strip()
            if txt:
                out.append((t, txt))
    seen, uniq = set(), []
    for t, p in out:                        # nested <sec> repeat paragraphs
        if p not in seen:
            seen.add(p)
            uniq.append((t, p))
    return uniq


def main():
    log = cfg.setup_logging("selection_phase6_lit")
    os.makedirs(OUT, exist_ok=True)
    cands = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(S, "results", "phase3",
                                                                       "candidates_annotated.csv")))}
    rows = []
    for pid in CASES:
        pmid = cands[pid]["pmid"]
        row = {"pdb_id": pid, "pmid": pmid, "citation": cands[pid]["citation"]}
        if not pmid:
            row["status"] = "no published paper linked (To be Published)"
            rows.append(row)
            log.info("%s: %s", pid, row["status"])
            continue
        q = urllib.parse.urlencode({"query": f"EXT_ID:{pmid} AND SRC:MED", "format": "json",
                                    "resultType": "core"})
        rec = json.loads(fetch(f"{EPMC}/search?{q}"))
        res = (rec.get("resultList") or {}).get("result") or []
        json.dump(res, open(os.path.join(OUT, f"{pid}_record.json"), "w"), indent=1)
        if not res:
            row["status"] = "no Europe PMC record"
            rows.append(row)
            continue
        r = res[0]
        row.update({"pmcid": r.get("pmcid", ""), "open_access": r.get("isOpenAccess", ""),
                    "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title", ""),
                    "year": r.get("pubYear", "")})
        excerpts = []
        if r.get("pmcid") and r.get("isOpenAccess") == "Y":
            try:
                xml = fetch(f"{EPMC}/{r['pmcid']}/fullTextXML")
                open(os.path.join(OUT, f"{pid}_fulltext.xml"), "w").write(xml)
                paras = paragraphs(xml)
                excerpts = [f"[{t}] {p}" for t, p in paras if TERMS.search(p)]
                row["status"] = f"open-access full text ({len(paras)} paragraphs, {len(excerpts)} relevant)"
            except Exception as e:
                row["status"] = f"full text download failed ({type(e).__name__}); abstract only"
        else:
            row["status"] = "full text not open access; abstract only"
        if not excerpts:
            abstract = re.sub(r"<[^>]+>", "", r.get("abstractText") or "")
            excerpts = [f"[abstract] {abstract}"] if abstract else []
        open(os.path.join(OUT, f"{pid}_excerpts.txt"), "w").write("\n\n".join(excerpts) + "\n")
        rows.append(row)
        log.info("%s: %s", pid, row["status"])
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(os.path.join(OUT, "literature_status.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
