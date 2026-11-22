"""Selection study, Phase 3b: annotate the frozen candidate list (without changing it).

Reads results/phase3/candidates_selected.csv (frozen; its SHA-256 is checked
against FROZEN_selection.txt) and writes candidates_annotated.csv with:
  why_useful         what the candidate tests, from its group (decision S6)
  literature_lead    curated reading of the abstract sentences found by the
                     keyword screen: "independent motion reported",
                     "rigidity reported", "linker/flexibility mentioned", or blank.
                     These notes were written by reading the sentences in
                     abstract_flexibility_sentences; they are leads to verify in
                     the full papers, not evidence on their own.
  notes              relations between candidates and caveats from metadata.

usage (from benchmarks/cyrange_comparison):
  python selection_study/scripts/phase3_annotate.py
"""
import csv
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.dont_write_bytecode = True
import bench_config as cfg  # noqa: E402

D = os.path.join(cfg.BENCH_DIR, "selection_study", "results", "phase3")

WHY = {
    "A": "two annotated domains: does RESUPER separate them where the reference merges them, "
         "and do they move relative to each other?",
    "B": "three or more annotated domains: the manuscript's strongest claim (reference merges "
         "multi-body arrangements)",
    "C": "negative control: one annotated domain with long unannotated (likely disordered) parts; "
         "RESUPER should not split the domain",
    "D": "known-weakness test: one domain plus a short separate helix (6EMR-like); does RESUPER "
         "miss a small rigid element?",
    "E": "over-splitting control: single-helix peptide (1A11/11TG/11TJ-like)",
}

# Curated from the abstract sentences (column abstract_flexibility_sentences).
LEADS = {
    "1B8T": "independent motion reported: two LIM domains are independent folding units with no apparent interaction",
    "1OQY": "independent motion reported: four structured domains connected by flexible linkers",
    "1Q8K": "independent motion reported: two domains mobile relative to each other",
    "2KDO": "linker/flexibility reported: conformational exchange in the N-terminal/central domain linker (NMR dynamics)",
    "2P01": "independent motion reported: three independent three-helix bundles joined by long flexible linkers",
    "6BYV": "independent motion reported: CR4 flexible relative to CR2-CR3, which interact",
    "1F5Y": "linker/flexibility reported: four-residue linker with low order parameters",
    "1WHE": "independent motion reported: Ca2+-dependent hinge between Gla and EGF domains",
    "2M8P": "independent motion reported: C-terminal domain orientation varies (RDC/SAXS)",
    "2RQ6": "linker/flexibility mentioned: flexibility of C-terminal domains examined by MD",
    "4A52": "rigidity reported: no preexisting conformational exchange in the apo enzyme (relaxation dispersion)",
    "5JS8": "state-dependent: rigid, closed arrangement with GTP analogues; mobility in other states",
    "2F3J": "linker/flexibility reported: transient N-terminal helix linked to the RRM by a flexible arm",
    "2N84": "two domains (FHA, PPIase) joined by a linker reported, although only one domain is annotated",
}

NOTES = {
    "2KDO": "same protein as 2L9N (SBDS): not independent",
    "2L9N": "same protein as 2KDO (SBDS): not independent",
    "6EMP": "sibling entry of 6EMR (same LEDGF/p75 study); 6EMR was excluded as RESUPER development data",
    "1F5Y": "homologous repeat protein to 6BYV (LDL/VLDL receptor ligand-binding repeats)",
    "6BYV": "homologous repeat protein to 1F5Y",
    "2M8P": "100 models from ensemble refinement against RDC and SAXS data, not a conventional NMR bundle",
    "1OQY": "only 12 models",
    "2N7L": "engineered chimera (troponin C/troponin I)",
    "7VH9": "engineered fusion construct",
    "7XEZ": "engineered fusion construct",
}


def main():
    path = os.path.join(D, "candidates_selected.csv")
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    frozen = open(os.path.join(D, "FROZEN_selection.txt")).read()
    if digest not in frozen:
        raise SystemExit("candidates_selected.csv no longer matches its frozen checksum")
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r["why_useful"] = WHY[r["group"]]
        r["literature_lead"] = LEADS.get(r["pdb_id"], "")
        note = NOTES.get(r["pdb_id"], "")
        if r["group"] == "D" and sum(int(r[f"n_dom_{s}"]) > 0 for s in ("CATH", "SCOP", "ECOD")) == 1:
            note = "; ".join(x for x in (note, "domain evidence from a single classification") if x)
        if not r["pmid"]:
            note = "; ".join(x for x in (note, "no PubMed-linked primary citation") if x)
        r["notes"] = note
    fields = list(rows[0].keys())
    with open(os.path.join(D, "candidates_annotated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"checksum OK ({digest[:16]}); wrote candidates_annotated.csv ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
