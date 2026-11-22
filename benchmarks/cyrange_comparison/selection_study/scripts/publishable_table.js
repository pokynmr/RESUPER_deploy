// Journal-style (three-rule) Word table of results/phase9/resuper_vs_cyrange_cases.json; converted to PDF.
// usage (from benchmarks/cyrange_comparison): NODE_PATH=$(npm root -g) node selection_study/scripts/publishable_table.js
const fs = require("fs"), path = require("path");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, BorderStyle, PageOrientation, AlignmentType } = require("docx");
const OUT = path.resolve(__dirname, "../results/phase9");
const D = JSON.parse(fs.readFileSync(path.join(OUT, "resuper_vs_cyrange_cases.json")));
const F = "Times New Roman";
const keep = ["PDB", "Protein", "Models / residues", "Domains RESUPER / CYRANGE", "RESUPER domains (RMSD)", "CYRANGE domains (RMSD)", "Evidence"];
const label = { "Models / residues": "Models / res.", "Domains RESUPER / CYRANGE": "Domains R / C" };
const W = [620, 1900, 760, 720, 2700, 2700, 5000];
const line = { style: BorderStyle.SINGLE, size: 6, color: "000000" }, none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const run = (t, o = {}) => new TextRun({ text: t, font: F, size: 15, ...o });
const cell = (t, k, o = {}) => new TableCell({ width: { size: W[k], type: WidthType.DXA }, columnSpan: o.span,
  borders: { top: o.top ? line : none, bottom: o.bottom ? line : none, left: none, right: none },
  margins: { top: 10, bottom: 10, left: 55, right: 55 },
  children: String(t).split(" | ").map((s, q, all) => new Paragraph({ children: [run(s, { bold: o.bold, italics: o.italic }),
    ...(o.ref && q === all.length - 1 ? [run(String(o.ref), { superScript: true })] : [])] })) });
const SEC = { A: "A. RESUPER separates domains that CYRANGE keeps in one range",
              B: "B. CYRANGE reports a region described as disordered or flexible as a separate domain; RESUPER does not (other domains identical)" };
const rows = [new TableRow({ tableHeader: true, children: keep.map((k, n) => cell(label[k] || k, n, { top: true, bottom: true, bold: true })) })];
let last = "";
D.rows.forEach((r, n) => {
  if (r.Category !== last) { last = r.Category;
    rows.push(new TableRow({ cantSplit: true, children: [cell(SEC[r.Category], 0, { span: keep.length, italic: true })] })); }
  rows.push(new TableRow({ cantSplit: true, children: keep.map((k, j) => cell(r[k], j, { bottom: n === D.rows.length - 1, bold: k === "PDB", ref: k === "Evidence" ? r.Ref : undefined })) }));
});
const nA = D.rows.filter((r) => r.Category === "A").length, nB = D.rows.length - nA;
const p = (runs, o = {}) => new Paragraph({ alignment: AlignmentType.JUSTIFIED, spacing: { after: 60 }, ...o, children: runs });
const t = (s, o = {}) => run(s, { size: 16, ...o });
const doc = new Document({ sections: [{ properties: { page: { size: { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE }, margin: { top: 560, right: 720, bottom: 500, left: 720 } } },
  children: [
    p([t("Table 1. ", { bold: true, size: 19 }), t(`NMR structures in which RESUPER's domain assignment is better supported than that of CYRANGE (${D.rows.length} entries: ${nA} in category A, ${nB} in category B).`, { size: 19 })], { spacing: { after: 100 } }),
    new Table({ width: { size: 14400, type: WidthType.DXA }, columnWidths: W, rows }),
    p([t("RESUPER: repository version of 2026-11-20, default parameters. CYRANGE: domains from each entry's wwPDB NMR validation report (identical to distributed CYRANGE 2.0 reports for 108/108 entries checked). RMSD: backbone (N, Cα, C′) RMSD to the mean structure, computed identically for both methods; a lower RMSD alone does not indicate a better assignment. Domains R / C: number of domains reported by RESUPER / CYRANGE. ", {}),
       t("Selection. ", { bold: true }), t("From 8,709 PDB NMR entries (solution NMR, ≥ 10 models, one protein chain, not used to develop RESUPER) compared with CYRANGE, entries in which the methods differed were judged by rules fixed before the scan; 164 entries favoured RESUPER on CATH/SCOP/ECOD annotation or ensemble rigidity and motion. Category A: every such entry in which RESUPER reports at least two domains, all of ≥ 20 residues and more than CYRANGE, with all RESUPER domains rigid (≤ 1.0 Å), no split inside one helix, and no contradiction in the original publication. Category B: entries of ≥ 60 residues in which every RESUPER domain matches a CYRANGE domain (Jaccard ≥ 0.9) and the remaining CYRANGE domains are non-rigid (> 1.0 Å, ≥ 8 residues), and the original publication explicitly describes that region as disordered or flexible; candidates were checked in order of decreasing RMSD. ")], { spacing: { before: 100, after: 60 } }),
    p([t("These entries were selected because they favour RESUPER and do not estimate how often RESUPER performs better. In 222 randomly drawn structures, RESUPER and CYRANGE gave identical domains for 136 and showed no significant difference in annotated domains recovered (162-structure confirmation sample: 32 vs 38, p = 0.125). In the scan, RESUPER reported more domains than CYRANGE, all of ≥ 20 residues, in only 17 entries; the evidence favoured CYRANGE in 10 of them (and one of the remaining 7, 2JT8, is excluded here because one RESUPER domain is not rigid).", { italics: true })]),
    ...[p([t("References", { bold: true })], { spacing: { before: 80, after: 20 } }),
      new Paragraph({ alignment: AlignmentType.JUSTIFIED, children: D.rows.flatMap((r, q) => { const [a, j, v, pg, y] = r._ref;
        return [t(`${q + 1}. ${a}, `), t(j, { italics: true }), t(" "), t(v, { bold: true }), t(`, ${pg} (${y})${q < D.rows.length - 1 ? "; " : "."}`)]; }) })],
  ] }] });
Packer.toBuffer(doc).then((b) => { fs.writeFileSync(path.join(OUT, "resuper_vs_cyrange_cases.docx"), b); console.log("docx ok"); });
