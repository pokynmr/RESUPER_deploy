"""Real CYRANGE domains from wwPDB NMR validation reports (decision S12).

The wwPDB validation pipeline runs CYRANGE (Kirchner & Guntert 2011) on every NMR
entry and stores, per residue of model 1, attribute cyrange_domain_id. We read
those, rebuild domain ranges, and compute each domain's backbone RMSD on OUR copy
of the structure with RESUPER's own functions (reference-mode atom convention).
The RMSD is therefore ours, not CYRANGE's printed value.

Validation: for PDB entries in the repository with a SAVED CYRANGE 2.0 report
(testdata/*.out, testdata/oracle_perturb/<ID>__ident.out), compare the wwPDB
domains with the saved report (residue numbers; exact match / same count).

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/wwpdb_cyrange.py validate
  python selection_study/scripts/wwpdb_cyrange.py cases       # the 22-row table
Outputs: selection_study/results/wwpdb_cyrange/
"""
import csv, gzip, os, re, sys, time, urllib.request
import xml.etree.ElementTree as ET
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts")); sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402
from report_parser import parse_report  # noqa: E402

S = os.path.join(cfg.BENCH_DIR, "selection_study")
OUT = os.path.join(S, "results", "wwpdb_cyrange")
XML = os.path.join(S, "data", "wwpdb_validation")
REPO = cfg.REPO_TESTDATA


def fetch(pid):
    os.makedirs(XML, exist_ok=True)
    p = os.path.join(XML, f"{pid.upper()}_validation.xml.gz")
    if not os.path.exists(p):
        l = pid.lower()
        url = f"https://files.rcsb.org/pub/pdb/validation_reports/{l[1:3]}/{l}/{l}_validation.xml.gz"
        for k in range(4):
            try:
                data = urllib.request.urlopen(url, timeout=60).read(); break
            except Exception as e:  # noqa: BLE001
                if k == 3 or "404" in str(e):
                    return None, f"download failed: {e}"
                time.sleep(2 ** k)
        open(p, "wb").write(data)
    return p, ""


def cyrange_domains(pid):
    """-> (status, n_domains_reported, {domain_id: set((chain, resnum))})"""
    p, err = fetch(pid)
    if not p:
        return err, None, {}
    # Fix 2026-11-20: use the FIRST model present, not literally model "1" (e.g. 1CRR numbers its
    # models 21-40; reading only model "1" silently gave 0 domains for 1CRR, 1EZC, 1EZD).
    bym, entry = {}, {}
    for _, el in ET.iterparse(gzip.open(p)):
        if el.tag == "Entry":
            entry = dict(el.attrib)
        elif el.tag == "ModelledSubgroup":
            a = el.attrib
            if a.get("cyrange_domain_id") and a.get("icode", " ").strip() == "" and (a.get("model") or "").isdigit():
                bym.setdefault(int(a["model"]), {}).setdefault(int(a["cyrange_domain_id"]), set()).add((a["chain"], int(a["resnum"])))
            el.clear()
    doms = bym[min(bym)] if bym else {}
    st = entry.get("cyrange_error", "no CYRANGE field")
    n = entry.get("cyrange_number_of_domains")
    return st, (int(n) if n not in (None, "") else None), dict(sorted(doms.items()))


def ranges(s):
    out, byc = [], {}
    for c, n in s:
        byc.setdefault(c, []).append(n)
    for c in sorted(byc):
        v = sorted(byc[c]); a = b = v[0]
        for x in v[1:]:
            if x == b + 1: b = x
            else: out.append(f"{a}..{b}"); a = b = x
        out.append(f"{a}..{b}")
    return ", ".join(out), sorted(byc)


class Rmsd:
    def __init__(self, extra_dirs=()):
        self.extra = list(extra_dirs)  # searched first (e.g. a fresh download folder)
        cfg.import_resuper()
        from resuper.core import _backbone_matrix
        from resuper.geometry import rmsd_to_mean_converged
        from resuper.pdb import parse_pdb
        self.bm, self.r, self.parse, self.c = _backbone_matrix, rmsd_to_mean_converged, parse_pdb, {}

    def __call__(self, pid, res):
        if pid not in self.c:
            cands = [os.path.join(d, f"{pid}.pdb") for d in self.extra] + [
                     os.path.join(S, "data", "raw_pdb", f"{pid}.pdb"),
                     os.path.join(S, "data", "expanded_pdb", f"{pid}.pdb"),
                     os.path.join(cfg.BENCH_DIR, "data", "fidelity_pdb", f"{pid}.pdb"),
                     os.path.join(REPO, f"{pid}.pdb")]
            found = [c for c in cands if os.path.exists(c)]
            if not found:
                raise FileNotFoundError(f"{pid}: no coordinate file in {cands}")
            s = self.parse(found[0])
            self.c[pid] = (s, {(r.chain, r.num): r.index for r in s.residues})
        s, ix = self.c[pid]
        idx = sorted(ix[r] for r in res if r in ix)
        if len(idx) < 3:
            return None
        ens, _ = self.bm(s, s.conformers, idx, exclude_cterm=True)
        return self.r(ens)


def jac(a, b):
    return len(a & b) / len(a | b) if a | b else None


def validate():
    """wwPDB CYRANGE vs saved CYRANGE 2.0 reports for the same PDB entries."""
    saved = {}
    for f in sorted(os.listdir(REPO)):
        m = re.fullmatch(r"(\d[A-Za-z0-9]{3})\.out", f)
        if m: saved[m.group(1).upper()] = os.path.join(REPO, f)
    for f in sorted(os.listdir(os.path.join(REPO, "oracle_perturb"))):
        m = re.fullmatch(r"(\d[A-Za-z0-9]{3})__ident\.out", f)
        if m: saved.setdefault(m.group(1).upper(), os.path.join(REPO, "oracle_perturb", f))
    rows = []
    for pid, path in saved.items():
        rep = parse_report(open(path).read())
        sv = [{n for _, n in d.residue_set("X")} for d in rep.domains]
        st, n, doms = cyrange_domains(pid)
        ww = [{n for _, n in s} for s in doms.values()]
        exact = sorted(map(sorted, sv)) == sorted(map(sorted, ww))
        pa, pb = set().union(*sv) if sv else set(), set().union(*ww) if ww else set()
        rows.append({"pdb_id": pid, "saved_report": os.path.relpath(path, REPO), "wwpdb_status": st,
                     "ndom_saved": len(sv), "ndom_wwpdb": len(ww), "same_count": len(sv) == len(ww),
                     "exact_same_domains": exact, "pooled_jaccard": round(jac(pa, pb), 3) if pa | pb else "",
                     "ranges_saved": " | ".join(d.text for d in rep.domains),
                     "ranges_wwpdb": " | ".join(ranges(s)[0] for s in doms.values())})
        print(pid, st, len(sv), len(ww), exact, flush=True)
    write(os.path.join(OUT, "validation_vs_saved_cyrange.csv"), rows)
    ok = [r for r in rows if r["wwpdb_status"] == "success"]
    multi = [r for r in ok if r["ndom_saved"] >= 2]
    s = (f"entries with saved CYRANGE report and a PDB ID: {len(rows)}; wwPDB CYRANGE success: {len(ok)}\n"
         f"same domain count: {sum(r['same_count'] for r in ok)}/{len(ok)}\n"
         f"exactly the same domains: {sum(r['exact_same_domains'] for r in ok)}/{len(ok)}\n"
         f"multi-domain (saved >= 2): same count {sum(r['same_count'] for r in multi)}/{len(multi)}, "
         f"exact {sum(r['exact_same_domains'] for r in multi)}/{len(multi)}\n")
    open(os.path.join(OUT, "validation_summary.txt"), "w").write(s); print(s)


def cases():
    calc = Rmsd()
    tab = list(csv.DictReader(open(os.path.join(S, "results", "phase7", "resuper_cases_table.csv"), encoding="utf-8-sig")))
    rows = []
    for t in tab:
        pid = t["PDB"]
        st, n, doms = cyrange_domains(pid)
        txt = []
        for s in doms.values():
            r, ch = ranges(s)
            v = calc(pid, s)
            txt.append(f"{r} ({v:.2f} Å)" if v is not None else r)
        rows.append({"PDB": pid, "Section": t["Source"], "Better (table)": t["Better"],
                     "Domains RESUPER": t["Domains RESUPER / comparator"].split(" / ")[0],
                     "Domains reference/saved": t["Domains RESUPER / comparator"].split(" / ")[1],
                     "Domains real CYRANGE (wwPDB)": len(doms) if st == "success" else "",
                     "wwPDB CYRANGE status": st,
                     "RESUPER domains (RMSD)": t["RESUPER domains (RMSD)"],
                     "Comparator domains (RMSD)": t["Comparator domains (RMSD)"],
                     "Real CYRANGE domains (RMSD computed here)": " | ".join(txt)})
        print(pid, st, len(doms), " | ".join(txt), flush=True)
    write(os.path.join(OUT, "cases_real_cyrange.csv"), rows)


def write(path, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    {"validate": validate, "cases": cases}[sys.argv[1]]()
