"""Publishable table (decision S15): RESUPER (repo_v2, default) vs CYRANGE (wwPDB validation reports).
Category A: RESUPER separates domains that CYRANGE merges. Category B: CYRANGE reports a non-rigid region
as a domain, RESUPER does not. Rules in docs/decisions_log.md (S15). Checks every rule from ranked.csv.
usage (from benchmarks/cyrange_comparison): python3 selection_study/scripts/publishable_table.py"""
import csv, json, os, re
B_ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
OUT = os.path.join(B_, "phase9"); os.makedirs(OUT, exist_ok=True)
rk = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(B_, "expanded", "ranked.csv")))}
pool = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(B_, "phase3", "candidate_pool.csv")))}
pp = {r["pdb_id"]: r for r in csv.DictReader(open(os.path.join(B_, "real_cyrange", "per_protein.csv")))}
def doms(s):
    out = []
    for d in s.split(" | "):
        if not d.strip(): continue
        res = set()
        for a, b in re.findall(r"(-?\d+)\.\.(-?\d+)", d): res |= set(range(int(a), int(b) + 1))
        m = re.search(r"\(([\d.]+) Å\)", d); out.append((res, float(m.group(1)) if m else None, d.strip()))
    return out
J = lambda a, b: len(a & b) / len(a | b)
A = {"2BW2": "ECOD: 2 vs 0 annotated domains recovered; CYRANGE range non-rigid; paper: two domains",
     "1BBO": "CATH: 2 vs 0; paper: RMSD reported per zinc finger",
     "2EW9": "CATH: 2 vs 0; paper: two ferredoxin folds, but they move as a unit (partial support)",
     "5MS9": "ECOD: 1 vs 0; CYRANGE range non-rigid; paper: domains joined by a flexible pivot",
     "6G81": "RESUPER-only domain rigid and moving; paper: Ni- and Zn-binding domains that reorient",
     "2V93": "RESUPER-only domain rigid and moving; paper: two domains, open/closed (PRE-refined ensemble)"}
Bc = {"9FLG": "disordered regions flank the helix-loop-helix (96-142)",
      "1HN6": "N- and C-terminal regions unstructured",
      "2MCD": "N-terminal PEST-like region unstructured",
      "7YWR": "helix fades after residue 50; 74-84 disordered",
      "1K9C": "residues 189-213 disordered",
      "2K8P": "N- and C-terminal arms unstructured, highly flexible",
      "2PXG": "N- and C-termini highly flexible",
      "2LWR": "N terminus disordered (A417-I429)",
      "2KTV": "C-terminal tail 414-437 flexible (15N relaxation)",
      "5NR5": "N- and C-terminal tails largely disordered"}
def name(p):
    n = (pool.get(p) or {}).get("protein") or (pp.get(p) or {}).get("protein") or ""
    n = n.capitalize() if n.isupper() else n
    return {"7YWR": n + " (Nsp8 N-terminal region)"}.get(p, n)
rows, problems = [], []
for cat, d in (("A", A), ("B", Bc)):
    for p, note in d.items():
        r = rk.get(p)
        if not r or r["qualifies"] != "True": problems.append(f"{p}: not qualifying after re-run"); continue
        R, C = doms(r["ranges_resuper"]), doms(r["ranges_cyrange"])
        if cat == "A":
            big = lambda D: sum(len(x[0]) >= 20 for x in D)
            ok = big(R) > big(C) and all(x[1] <= 1.0 for x in R) and r["helix_split"] != "True"
        else:
            used = {max(range(len(C)), key=lambda j: J(x[0], C[j][0])) for x in R}
            match = all(max(J(x[0], c[0]) for c in C) >= 0.9 for x in R)
            extra = [C[j] for j in range(len(C)) if j not in used]
            ok = match and extra and all(e[1] > 1.0 and len(e[0]) >= 8 for e in extra) and int(r["n_residues"]) >= 60
        if not ok: problems.append(f"{p}: rule check failed"); continue
        if cat == "B":
            note = f"Extra CYRANGE domain{'s' if len(extra) > 1 else ''} {', '.join(f'{e[1]:.2f}' for e in extra)} Å; paper: {note}"
        rows.append({"Category": cat, "PDB": p, "Protein": name(p), "Models / residues": f"{r['n_models']} / {r['n_residues']}",
                     "Domains RESUPER / CYRANGE": f"{r['ndom_resuper']} / {r['ndom_cyrange']}",
                     "RESUPER domains (RMSD)": r["ranges_resuper"].replace("..", "–"),
                     "CYRANGE domains (RMSD)": r["ranges_cyrange"].replace("..", "–"),
                     "Evidence": re.sub(r"(\w)-(\w?\d)", r"\1–\2", note).replace("15N", "¹⁵N"), "Ref": len(rows) + 1})
JN = {"J Biol Chem": "J. Biol. Chem.", "Proc Natl Acad Sci U S A": "Proc. Natl. Acad. Sci. USA", "J Biol Inorg Chem": "J. Biol. Inorg. Chem.",
      "J Mol Biol": "J. Mol. Biol.", "Nucleic Acids Res": "Nucleic Acids Res.", "Nat Commun": "Nat. Commun.", "Sci Adv": "Sci. Adv."}
refs = json.load(open(os.path.join(OUT, "refs", "refs.json")))
def cite(m):
    au = [a.strip() for a in m["authors"].rstrip(".").split(",")]
    def fmt(a):
        parts = a.split(); last, ini = " ".join(parts[:-1]), parts[-1]
        return " ".join(c + "." for c in ini) + " " + last
    names = fmt(au[0]) + " et al." if len(au) > 3 else ", ".join(fmt(a) for a in au[:-1]) + (" and " if len(au) > 1 else "") + fmt(au[-1])
    pg = (m["pages"] or "").split("-")[0]
    return names, JN.get(m["journal"], m["journal"]), m["volume"], pg, m["year"]
for r in rows:
    n, j, v, pg, y = cite(refs[r["PDB"]]); r["Reference"] = f"{n}, {j} {v}, {pg} ({y})"; r["_ref"] = [n, j, v, pg, y]
print("rows", len(rows)); print("problems", problems)
F = [k for k in rows[0] if not k.startswith("_")]
with open(os.path.join(OUT, "resuper_vs_cyrange_cases.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=F, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
json.dump({"fields": F, "rows": rows, "problems": problems}, open(os.path.join(OUT, "resuper_vs_cyrange_cases.json"), "w"), indent=1, ensure_ascii=False)
