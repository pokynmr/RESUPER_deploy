"""Step 2b: inspect every structure and complete the dataset manifest.

For each PDB ID in data/dataset_manifest_DRAFT.csv this script finds the
structure file (data/raw_pdb/<ID>.pdb; for IDs not downloaded yet it falls back
to the RESUPER repository's testdata/ copy, read only) and records:

  from the raw file      models, chains, residue numbering and jumps, insertion
                         codes, alternate locations, non-standard / HETATM
                         residues, whether all models contain the same atoms,
                         experimental method, molecule name
  from RESUPER's reader  models and residues RESUPER will actually use, and
                         residues it drops (missing N, CA or C in some model)
  from RCSB metadata     release / last revision date and status (if fetched)

The difference between the two views matters: RESUPER silently skips
non-standard residues, merges residues that differ only by insertion code, and
drops residues with incomplete backbones. Any such case is listed in `issues`.

If both an RCSB download and a repository copy exist (1A67, 1ADR), their
coordinates are compared. The saved CYRANGE reports were made on the
repository copies, so we need to know whether the two are identical.

Outputs (the draft manifest is never modified):
  data/structure_inventory.csv   one row per structure, all measurements
  data/dataset_manifest.csv      the draft manifest with TBD fields filled in
  logs/inspect_<time>.log

usage (from benchmarks/cyrange_comparison):
  python scripts/inspect_structures.py
"""
import csv
import hashlib
import json
import os
from collections import OrderedDict

import bench_config as cfg

STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
IGNORE_HET = {"HOH", "DOD", "WAT"}


def scan_raw(path):
    """Read the PDB file line by line, independent of RESUPER's parser."""
    info = {"title": [], "molecules": [], "expdta": "", "revdat_latest": ""}
    models = []                 # per model: set of (chain, resnum, icode, resname, atom)
    residues = OrderedDict()    # model-1 residues in file order -> record type
    altlocs, het = set(), set()
    cur = None
    coord_hash = hashlib.sha256()
    with open(path, errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec.startswith("TITLE"):
                info["title"].append(line[10:80].strip())
            elif rec.startswith("COMPND") and "MOLECULE:" in line:
                info["molecules"].append(line.split("MOLECULE:", 1)[1].strip().rstrip(";"))
            elif rec.startswith("EXPDTA"):
                info["expdta"] = (info["expdta"] + " " + line[10:80].strip()).strip()
            elif rec.startswith("REVDAT") and not info["revdat_latest"]:
                info["revdat_latest"] = line[13:22].strip()  # first REVDAT = newest
            elif rec.startswith("MODEL"):
                cur = set()
                models.append(cur)
            elif rec in ("ATOM  ", "HETATM"):
                if cur is None:                 # file without MODEL records
                    cur = set()
                    models.append(cur)
                resname = line[17:20].strip()
                if rec == "HETATM" and resname in IGNORE_HET:
                    continue
                key = (line[21].strip() or "A", line[22:26].strip(), line[26].strip(),
                       resname)
                cur.add(key + (line[12:16].strip(),))
                if line[16].strip():
                    altlocs.add(line[16])
                if rec == "HETATM":
                    het.add(resname)
                if len(models) == 1 and key not in residues:
                    residues[key] = rec
                # model index + atom identity + coordinates (cols 13-54)
                coord_hash.update(f"{len(models)}|{line[12:54]}\n".encode())
    info["title"] = " ".join(info["title"])
    info["models"] = models
    info["residues"] = residues
    info["altlocs"] = altlocs
    info["het"] = het
    info["coord_sha256"] = coord_hash.hexdigest()
    return info


def numbering_summary(residues):
    """Per chain: 'A:9..116' plus any jumps in numbering (e.g. 'A:47->49')."""
    chains = OrderedDict()
    for (ch, num, ic, resname), rec in residues.items():
        # Only protein residues count: ions, caps and ligands (HETATM or
        # non-standard names) have their own numbers and are not chain gaps.
        if rec != "ATOM  " or resname not in STANDARD:
            continue
        try:
            chains.setdefault(ch, []).append((int(num), ic))
        except ValueError:
            pass
    spans, jumps = [], []
    for ch, nums in chains.items():
        spans.append(f"{ch}:{nums[0][0]}..{nums[-1][0]}")
        for (a, _), (b, ic) in zip(nums, nums[1:]):
            if b != a + 1 and not (b == a and ic):
                jumps.append(f"{ch}:{a}->{b}")
    return chains, spans, jumps


def load_metadata(pdb_id):
    path = os.path.join(cfg.META_DIR, f"{pdb_id}.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        m = json.load(fh)
    acc = m.get("rcsb_accession_info", {})
    return {
        "rcsb_title": m.get("struct", {}).get("title", ""),
        "rcsb_method": ";".join(e.get("method", "") for e in m.get("exptl", [])),
        "rcsb_initial_release": acc.get("initial_release_date", "")[:10],
        "rcsb_revision_date": acc.get("revision_date", "")[:10],
        "rcsb_status": m.get("pdbx_database_status", {}).get("status_code", ""),
        "rcsb_model_count": m.get("rcsb_entry_info", {}).get("deposited_model_count", ""),
    }


def find_files(pdb_id):
    """Return [(label, path)] for the RCSB download and/or the repository copy."""
    out = []
    raw = os.path.join(cfg.RAW_PDB_DIR, f"{pdb_id}.pdb")
    if os.path.exists(raw):
        out.append(("rcsb_download", raw))
    if os.path.isdir(cfg.REPO_TESTDATA):
        for f in os.listdir(cfg.REPO_TESTDATA):
            if f.lower() == f"{pdb_id.lower()}.pdb":
                out.append(("repo_testdata", os.path.join(cfg.REPO_TESTDATA, f)))
    return out


def inspect(pdb_id, label, path, parse_pdb, log):
    row = {"pdb_id": pdb_id, "file_role": label, "path": path}
    raw = scan_raw(path)
    st = parse_pdb(path)
    chains, spans, jumps = numbering_summary(raw["residues"])
    icodes = [f"{c}:{n}{i}" for (c, n, i, _) in raw["residues"] if i]
    nonstd = sorted({r for (_, _, _, r), rec in raw["residues"].items()
                     if rec == "ATOM  " and r not in STANDARD})
    ref_atoms = raw["models"][0] if raw["models"] else set()
    n_diff = sum(1 for m in raw["models"][1:] if m != ref_atoms)
    backbone = {"N", "CA", "C"}
    ref_bb = {a for a in ref_atoms if a[4] in backbone}
    n_diff_bb = sum(1 for m in raw["models"][1:]
                    if {a for a in m if a[4] in backbone} != ref_bb)
    n_raw_res = sum(1 for (_, _, _, r), rec in raw["residues"].items()
                    if rec == "ATOM  " and r in STANDARD)
    row.update({
        "n_models_file": len(raw["models"]),
        "n_models_resuper": st.n_models,
        "n_residues_file_standard": n_raw_res,
        "n_residues_resuper": len(st.residues),
        "chains": " ".join(f"{c}({len(v)})" for c, v in chains.items()),
        "numbering": " ".join(spans),
        "numbering_jumps": " ".join(jumps),
        "insertion_codes": " ".join(icodes),
        "altlocs": "".join(sorted(raw["altlocs"])),
        "nonstandard_in_ATOM": " ".join(nonstd),
        "hetatm_residues": " ".join(sorted(raw["het"])),
        "models_with_different_atoms": n_diff,
        "models_with_different_backbone": n_diff_bb,
        "resuper_dropped": " ".join(f"{c}{n}" for c, n in st.dropped),
        "expdta": raw["expdta"],
        "molecules": " | ".join(dict.fromkeys(raw["molecules"])),
        "title": raw["title"],
        "revdat_latest": raw["revdat_latest"],
        "coord_sha256": raw["coord_sha256"],
    })
    row.update(load_metadata(pdb_id))

    issues = []
    if st.n_models < 2:
        issues.append("fewer than 2 models: not an ensemble, RESUPER cannot run")
    if "NMR" not in raw["expdta"].upper():
        issues.append(f"experimental method is '{raw['expdta']}', not NMR")
    if len(chains) > 1:
        issues.append(f"{len(chains)} chains: ranges are chain-prefixed; check domains across chains")
    if icodes:
        issues.append("insertion codes present: RESUPER ignores them, so residues "
                      "sharing a number are MERGED")
    if raw["altlocs"]:
        issues.append("alternate locations present: RESUPER keeps the first")
    if nonstd or raw["het"]:
        issues.append("non-standard/HETATM residues skipped by RESUPER (may create gaps)")
    if jumps:
        issues.append("residue numbering is not consecutive (counted as gaps by RESUPER)")
    if n_diff_bb:
        issues.append(f"{n_diff_bb} model(s) differ from model 1 in BACKBONE atoms")
    elif n_diff:
        issues.append(f"{n_diff} model(s) differ from model 1 in non-backbone atoms "
                      "only (e.g. hydrogens; minor)")
    if st.dropped:
        issues.append(f"RESUPER drops {len(st.dropped)} residue(s) lacking N/CA/C")
    if n_raw_res != len(st.residues) + len(st.dropped) and not icodes:
        issues.append("residue count differs between raw file and RESUPER parser")
    row["issues"] = "; ".join(issues) if issues else "none found"
    row["eligible"] = "yes" if st.n_models >= 2 and len(st.residues) > 0 else "no"
    log.info("%s [%s]: %d models, %d residues, chains %s; %s", pdb_id, label,
             st.n_models, len(st.residues), row["chains"], row["issues"])
    return row


def main():
    log = cfg.setup_logging("inspect")
    parse_pdb = cfg.import_resuper().pdb.parse_pdb
    with open(cfg.MANIFEST_DRAFT) as fh:
        manifest = list(csv.DictReader(fh))

    inventory, by_id = [], {}
    for m in manifest:
        pdb_id = m["pdb_id"].upper()
        files = find_files(pdb_id)
        if not files:
            log.warning("%s: no structure file yet", pdb_id)
            continue
        rows = []
        for label, path in files:
            try:
                rows.append(inspect(pdb_id, label, path, parse_pdb, log))
            except Exception as e:  # record and continue
                log.error("%s [%s]: inspection failed: %s: %s", pdb_id, label,
                          type(e).__name__, e)
                rows.append({"pdb_id": pdb_id, "file_role": label, "path": path,
                             "issues": f"inspection failed: {type(e).__name__}: {e}",
                             "eligible": "no"})
        if len(rows) == 2 and rows[0].get("coord_sha256") and rows[1].get("coord_sha256"):
            same = rows[0]["coord_sha256"] == rows[1]["coord_sha256"]
            note = ("RCSB and repository copies have identical coordinates" if same else
                    "RCSB and repository copies DIFFER: saved CYRANGE report was made "
                    "on the repository copy")
            for r in rows:
                r["copy_comparison"] = note
            log.info("%s: %s", pdb_id, note)
        inventory.extend(rows)
        by_id[pdb_id] = rows[0]      # prefer the RCSB download for the manifest

    fields = list(OrderedDict.fromkeys(k for r in inventory for k in r))
    with open(cfg.INVENTORY, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(inventory)

    # Manifest: copy the draft and fill in what was measured.
    for m in manifest:
        r = by_id.get(m["pdb_id"].upper())
        if m["paper_set"] == "A":
            m["ensemble_matches_paper"] = ("NO - deposited-entry stand-in for the "
                                           "paper's recalculated bundles (decision D2)")
        if not r or "n_models_resuper" not in r:
            m["eligible"] = "pending (no structure file)" if not r else "no"
            continue
        m["n_models"] = r["n_models_resuper"]
        m["n_residues"] = r["n_residues_resuper"]
        m["chains"] = r["chains"]
        m["residue_numbering"] = r["numbering"] + (
            f" (jumps: {r['numbering_jumps']})" if r["numbering_jumps"] else "")
        m["local_status"] = f"present ({r['file_role']})"
        if r["file_role"] == "rcsb_download":
            m["structure_source"] = "RCSB deposited entry, data/raw_pdb/" + m["pdb_id"].upper() + ".pdb"
        if m["protein_name"].startswith("TBD") and r.get("molecules"):
            m["protein_name"] = r["molecules"]
            m["name_source"] = "PDB file COMPND"
        m["issues"] = "; ".join(x for x in (r["issues"], r.get("copy_comparison", ""),
                                            m["issues"]) if x and x != "none found") \
            or "none found"
        m["eligible"] = r["eligible"]
    fields = list(manifest[0].keys()) + (["eligible"] if "eligible" not in manifest[0] else [])
    with open(cfg.MANIFEST, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(manifest)
    n_present = sum(1 for m in manifest if m.get("eligible") in ("yes", "no"))
    log.info("inspected %d of %d manifest entries; wrote %s and %s", n_present,
             len(manifest), cfg.INVENTORY, cfg.MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
