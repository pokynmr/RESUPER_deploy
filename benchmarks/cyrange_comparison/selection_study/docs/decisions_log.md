# Selection study: decisions log

## S1: Which repository structures to investigate (2026-11-11, user + rule)
The user asked for "the structures in the repository that show RESUPER
identifying regions CYRANGE does not". Rather than take the manuscript's
examples, all 110 real ensembles with a saved real CYRANGE 2.0 report (rebuilt in
the earlier benchmark) were screened. Kept: any ensemble where RESUPER default
and CYRANGE differ in domain count, have an unmatched or weakly matched
(J < 0.5) domain, or a split/merge, in EITHER direction. Result: 1A11, 11TG,
11TJ, 1DX7, 6EMR, 1BFW. Cases going against RESUPER (6EMR) are kept on purpose.

## S2: Size-matched baseline for relative rotation (2026-11-11)
Rotation angles of short segments are noisier than of long ones, so the
earlier 12° baseline (protein halves) is not fair for 8–13-residue segments.
For each case, the baseline is built from segment pairs with the same lengths
and spacing as RESUPER's two domains, taken from inside the 29 proteins that
both methods call one rigid domain in the earlier benchmark. The observed
median angle is reported as a percentile of the baseline per-pair medians.

## S3: Evidence categories
Every conclusion is labelled: computational difference / supported by the
ensemble geometry / supported by independent (literature or experimental)
evidence / hypothesis. "Supported by the ensemble" does not mean "real
dynamics": NMR ensemble spread can reflect missing restraints.

## S4: Independent domain evidence (2026-11-11)
CATH, SCOP and ECOD domain annotations from the RCSB Data API
(polymer_entity_instance) are used as independent evidence of multiple domains.
Segments are grouped by their *_DOMAIN_ID; sequential positions are converted
to author numbering with auth_to_entity_poly_seq_mapping. Raw responses are kept
in results/phase2/annotations/raw/. Classification domains are folding units,
not proof of independent motion (1DX7).

## S5: RESUPER decision traces (2026-11-11)
scripts/phase2_split_trace.py replays RESUPER's Stage A walk with RESUPER's own
functions and must reproduce `_rmsd_split_groups` exactly before its log is used
(it did for 1A11, 11TG, 11TJ, 1DX7, 6EMR).

## S6: Phase 3 candidate selection rules (2026-11-11; fixed BEFORE the search was run)
User decisions: 10-15 proteins per group; group E included. Chosen: 12 per group.
Source: RCSB Search API (pool) + RCSB GraphQL Data API (annotations, metadata);
no coordinates downloaded and no method run during selection.
Common filters: experimental method SOLUTION NMR; >= 10 deposited models; exactly
one protein entity and one polymer instance (monomer), no nucleic acid.
Excluded before sampling: every PDB ID that appears in the RESUPER repository's
development material (testdata/ file names incl. oracle_* folders,
testdata/pacsy_ids.txt, scripts/bench/proteins.txt; these were used to build or
tune RESUPER) and the 37 CYRANGE-paper proteins of the earlier benchmark.
Domain counts: per source (CATH, SCOP, ECOD), segments grouped by *_DOMAIN_ID.
Positions are sequential (label_seq_id). L = sample sequence length.
  B  >= 2 sources annotate; every annotating source has >= 2 domains; median count >= 3; L 120-500
  A  same, but median count == 2; L 80-500
  D  >= 1 source; every source 1 domain; a helix (HELIX_P) of 8-25 residues lying
     entirely outside every source's domain, >= 3 residues from it; L 60-400
  C  >= 2 sources; every source 1 domain; residues outside the union of all
     domains >= max(20, 0.2 L); L 60-300
  E  L 25-50; longest helix covers >= 70% of L (classification not required)
  Exclusive assignment in the order B, A, D, C, E.
Sampling: 12 per group, uniformly at random with numpy seed 20261111 from the
eligible pool (fewer if the pool is smaller). The selected list is written to
results/phase3/candidates_selected.csv and frozen (SHA-256 recorded) before Phase 4.
Literature: for each selected entry with a PubMed ID, the abstract (Europe PMC)
is scanned for flexibility-related terms; matching sentences are recorded as
leads, not as evidence.

## S7: Phase 4 screening rules (2026-11-11; fixed BEFORE downloading or running)
Input: the 60 frozen candidates (checksum verified first). Deposited PDB-format
files from RCSB, stored unmodified with SHA-256. Both methods are run on every
file RESUPER can read, with default parameters, via the unmodified command line;
a second full run must reproduce every report.
Each structure is classified:
  SOFTWARE COMPATIBILITY problem (the structure may be fine; the tool cannot
  handle it as deposited):
    - no PDB-format file available (RESUPER reads PDB format only);
    - RESUPER exits with an error or exceeds 1 h;
    - insertion codes (RESUPER merges residues sharing a number);
    - non-standard/modified residues inside the chain, which RESUPER skips
      (<= 20% of residues: kept, flagged as "artificial gaps possible").
  SCIENTIFICALLY UNSUITABLE (excluded from Phase 5 statistics):
    - > 20% of polymer residues non-standard (e.g. D-amino-acid peptides);
    - RESUPER reads fewer than 20 residues;
    - models differ in backbone atoms (inconsistent ensemble);
    - more than one protein chain read.
  Otherwise SUITABLE (minor caveats such as missing hydrogens are recorded).
Results of excluded or failed structures are still saved, not deleted.

## S8: Phase 5 comparison rules (2026-11-12; fixed BEFORE the comparison was run)
A = RESUPER default, B = reference mode. Residues are (chain, author number).
1. A vs B: the earlier benchmark's procedure (Jaccard; Hungarian one-to-one
   matching; identical / strong >= 0.8 / partial >= 0.5 / weak; merge = a domain
   holding >= 50% of each of >= 2 domains of the other method; coverage counts a
   residue once). Protein label: identical / minor / substantial (as before).
2. Independent reference: per entry, the classification (CATH, ECOD, SCOP) whose
   domains cover the most residues; ties -> CATH > ECOD > SCOP. Annotated domains
   are converted to author numbering. Group E has no annotation.
   For each method:
   - mixed domain: >= 3 residues AND >= 10% of its residues from each of >= 2
     annotated domains;
   - recovered annotated domain: one method domain covers >= 50% of it and that
     domain is not mixed;
   - split annotated domain: >= 2 method domains each cover >= 20% of it;
   - missed annotated domain: no method domain covers >= 20% of it.
3. Relative motion: for every pair of domains reported by either method, the
   per-model rotation of domain 2 after superposing domain 1 (median), compared
   with a size-matched baseline from inside the 29 single-domain proteins of the
   earlier benchmark. Baseline segment lengths are capped at 40 residues and the
   gap at 20 (conservative: shorter segments are noisier). Fewer than 20
   baseline pairs -> "not calculable". "Supported by motion" = observed median
   above the baseline's 95th percentile.
4. Pattern labels per protein (from A vs B): P1 RESUPER more domains; P2 reference
   more domains; P3 same count, substantially different boundaries; identical;
   minor. A RESUPER domain with no reference counterpart (J = 0) and < 40
   residues counts as "smaller region found by RESUPER"; the reverse as "smaller
   region found by reference".
5. Evidence level of each difference: computational only / supported by motion /
   supported by annotation / supported by both.
6. Aggregates exclude 1WHE (severe modified-residue gaps) and are reported per
   group; non-independent pairs (2KDO/2L9N, 1F5Y/6BYV) are flagged.

## S8a: Amendment to S8.3, baseline caps (2026-11-12; after the first run, before interpretation)
Problem: with caps of 40 residues per segment and a gap of 20, each baseline pair
needs ~100 contiguous residues inside one rigid domain. 41 of 61 domain pairs
got fewer than 20 baseline pairs and were "not calculable", including pairs
rotating by 84-127 degrees (1OQY). The rule made the test uncomputable for most
proteins, independent of which method reported the domains.
Change: caps of 20 residues per segment and a gap of 10. Shorter segments are
noisier, so the baseline is HIGHER, making "supported by motion" harder to reach
(conservative). Minimum of 20 baseline pairs unchanged.
Also fixed (code bug, not a rule change): "supported_by_motion" was stored as a
NumPy boolean, so the summary's `is True` checks counted 0 supported pairs.

## S8b: Amendment to S8.2, inconsistent annotations (2026-11-12; before interpretation)
2LY2: CATH lists a domain over the whole chain (1-180) and a second inside it
(100-180); the chosen source's domains overlapped by 81 residues. Rule added:
a source whose domains overlap each other is treated as inconsistent and the
next source in the tie order (CATH > ECOD > SCOP) with non-overlapping domains is
used (2LY2 -> ECOD: 1-88, 89-180). Only 2LY2 of the 60 entries is affected.

## S9: Phase 6 scope (2026-11-12, user)
Investigated in depth: 1OQY, 2P01, 2KUI, 2BW2, 2RPC (RESUPER separates more),
2F3J, 7A2D (reference finds small elements RESUPER misses), 2LQ4 (RESUPER split
without motion), 2M70 (RESUPER control split that moves), 1B8T (reference
overlapping small pieces). The choice was made by judgement (proposed by Claude,
accepted by the user), not by rule. The other substantial disagreements (1Q8K,
2LY2, 6H5H, 5JS8, 8H0U, 8RPQ) and 1WHE are not investigated in depth; this is
stated as a limitation in the final report. The 7 minor differences are covered
by the known last-residue atom-counting convention.

## S10: Phase 7 ranking rubric (2026-11-12; fixed before ranking)
Purpose: prioritise structures for further investigation, not a statistical test.
Pool: the 10 Phase 6 proteins plus the Phase 1 repository cases with real
CYRANGE 2.0 output (1A11, 11TG, 11TJ, 1DX7, 6EMR).
Score for "RESUPER's division is the more meaningful one":
  +3 full-text literature supports RESUPER's division quantitatively (e.g.
     published domain boundaries, measured disorder of what RESUPER omits)
  +2 literature supports it qualitatively (abstract or title)
  +1 independent annotation (CATH/ECOD/SCOP) supports RESUPER over the comparator
  +1 ensemble geometry supports it (relative motion beyond the size-matched
     baseline AND internally rigid domains, or omitted region not rigid)
  +1 the comparator is the REAL CYRANGE 2.0 output (not reference mode)
  -1 per major caveat (ensemble spread partly from modelling; peptide with
     possible under-restraint; < 10 models; both methods share the main error)
  -2 literature supports the comparator's division instead
Ties: broken by how much of the protein the comparator leaves out or merges.
Structures with score <= 0, or where the comparator is better, are listed as
counter-examples or open cases, not ranked as promising.

## S11: Confirmation sample, pre-registration (2026-11-13; fixed BEFORE drawing)
User decisions: group A 40-50 (chosen: 45); group B: complete the census; also
expand the control groups (C, D, E).
Sample: drawn from results/phase3/candidate_pool.csv (the Phase 3 pool, groups by
rule S6), excluding the 60 Phase 3 candidates. numpy seed 20261113.
  A: 45 at random | B: all remaining (3) | D: all remaining (28)
  C: 45 at random | E: 45 at random            total 166
The list is frozen (SHA-256) before any structure is downloaded.
Screening (S7), runs (default parameters, repeat check) and comparison
(S8 with S8a, S8b) are exactly as before. Structures with >= 5 in-chain
non-standard residues (severe) or classed scientifically unsuitable are
excluded from the tests below.
Hypotheses and tests (two-sided, alpha = 0.05), on the confirmation sample:
  H1 (groups A+B): per protein, d = annotated domains recovered by RESUPER minus
     those recovered by reference mode. Wilcoxon signed-rank test on d.
  H2 (groups A+B): same for annotated domains missed (reference minus RESUPER).
  H3 (groups C+E, over-splitting): per protein, does each method report >= 2
     domains? Exact McNemar test (binomial test on discordant proteins).
  H4 (group D, short elements): does each method report a domain holding >= 50%
     of the residues of the short helix that defined the group? Exact McNemar.
Descriptive (no test): two-part rigid-domain check for every domain reported by
only one method: internal RMSD <= 1.0 A AND relative rotation above the
size-matched baseline's 95th percentile (motion alone is not enough: Phase 6).
Reporting: original 60 and confirmation 166 reported separately, then combined.
Comparator remains reference mode: conclusions concern RESUPER vs reference mode.

## S12: real CYRANGE from wwPDB validation reports; papers for the 10 showcase cases
- Source: the wwPDB NMR validation XML (files.rcsb.org) stores CYRANGE's domain
  per residue (`cyrange_domain_id`, model 1). Script `scripts/wwpdb_cyrange.py`.
- Check before use: for all 108 PDB entries with a saved CYRANGE 2.0 report in the
  repository (testdata), wwPDB domains are exactly identical (108/108; 9/9 multi-domain).
  2 entries had no wwPDB CYRANGE result (1ACP, 1BFW) and were not used.
- 22 table cases: real CYRANGE equals RESUPER for 9 (1OQY within 4 residues); of the
  17 judged RESUPER-better against reference mode, 10 are not differences from real CYRANGE.
- Papers (Europe PMC, `scripts/s12_literature.py`) describe 2KDR and 9G7T as a single helix,
  first read as contradicting RESUPER's split (revised in S16: not a contradiction) and favour CYRANGE's split for 7YYI (GxxxG hinge). 7A2D: CYRANGE separates
  the two BON domains; marked unclear.
- Outputs: results/wwpdb_cyrange/, results/phase7/resuper_cases_table_v2.{csv,docx,pdf}.
  Earlier table (v1) kept unchanged.

## S13: RESUPER vs real CYRANGE on all 222 proteins
- Script `scripts/real_cyrange_all.py`; same comparison functions and H1-H4 rules as S11,
  real CYRANGE (wwPDB, S12) as comparator. Post hoc (H1-H4 were registered against reference mode).
- Result: no significant difference on any test; RESUPER identical to CYRANGE for 136/222.
  Write-up: docs/real_cyrange_results.md.

## S14: expanded scan of the whole eligible pool against real CYRANGE (written before the scan)
Purpose (user request): find 20 showcase proteins where RESUPER finds more domains or
works better than REAL CYRANGE. These are chosen BECAUSE they favour RESUPER: a
showcase, not a benchmark; the unbiased comparison is S13 (222 proteins).
- Pool: Phase 3 pool (9,942 entries: solution NMR, >= 10 models, one protein entity and
  instance, no RESUPER development entries), minus the 226 already drawn; all of them,
  grouped or not. RESUPER default parameters only. Real CYRANGE from wwPDB (S12).
- Scan: per entry, wwPDB CYRANGE first (no result: skipped); then download, S7
  screening, RESUPER run (timeout 600 s). Structures are deleted unless RESUPER and
  CYRANGE differ (domain count or "substantial difference"), to save disk.
- Exclusions: SCIENTIFICALLY UNSUITABLE (S7), >= 5 in-chain non-standard residues.
- Evidence (only for entries that differ), same functions as Phase 5 / S11:
  E1 annotation (CATH/SCOP/ECOD, S8b rule): RESUPER recovers more annotated domains than
     CYRANGE, or as many with fewer mixed domains.
  E2 a CYRANGE domain is not internally rigid (RMSD > 1.0 A) while all RESUPER domains are.
  E3 a RESUPER-only domain passes the two-part check (rigid <= 1.0 A AND moves relative
     to the rest).
  Disqualifying: X1 annotation favours CYRANGE; X2 a CYRANGE-only domain passes the
  two-part check; X3 a RESUPER-only domain is not internally rigid.
- Ranking: tier 1 = E1; tier 2 = E2; tier 3 = E3 only; tier 4 = E3 only where RESUPER's
  split lies inside one deposited HELIX record (helix split; see S16 for why this is
  not evidence against the split). Within a tier: larger annotation gain, then larger CYRANGE RMSD
  excess, then motion percentile.
- Literature: abstract (and open full text) of the top 20 checked; an entry whose paper
  explicitly contradicts RESUPER's division is replaced by the next-ranked one and listed.

### S14 outcome (2026-11-17)
- Expanded scan: 9,716 entries; 8,709 compared with real CYRANGE; 693 excluded (S7 screening); 303 no wwPDB CYRANGE result; 11 not downloadable; 862 disagreements kept.
- Evidence step (expanded_evidence.py): 908 candidates (862 + earlier differing entries); 165 qualify (tier 1: 37, 2: 93, 3: 24, 4: 11).
- Literature (s14_literature.py, top 30): 1K8V replaced (paper: residues up to 13 random coil; RESUPER domain 9..32 includes them, 2.54 Å). 1WIL: citation search returned an unrelated article, kept with "no relevant paper". 2EW9 noted as partly supported (two ferredoxin folds, but tumble as a unit).
- Table: results/phase8/showcase20_table.{csv,json,docx,pdf} (scripts showcase20_table.py, showcase20_pdf.js). All 20 are tier 1; ranks 7-21 are mostly coverage differences in single-domain peptides, stated in the table note. Labelled as a selection favouring RESUPER, with the S13 unbiased result alongside.

### S14b multi-domain showcase (2026-11-18, user request: replace ranks 7-21 with stronger multi-domain cases)
- Criterion fixed before paper check: qualifying under S14, not a helix split, RESUPER reports >=2 domains each >=20 residues. Ranks 1-6 of showcase20 kept as the user asked.
- Result: 10 entries meet it (2BW2, 1BBO, 2EW9, 6BYV, 5MS9, 1B8T, 2KUO, 2JT8, 6G81, 2V93) + 1CRR kept = 11 rows. Not padded to 20.
- Across all 908 disagreements, RESUPER reports more >=20-residue domains than CYRANGE in only 17 entries: 7 qualify, 10 excluded (evidence favours CYRANGE).
- Papers: 1B8T, 6G81, 2V93 support RESUPER; 2KUO no relevant paper found; 2JT8 weaker (N-lobe split not supported, one range 1.87 Å). 2KLO, 2YRQ near-misses (extra domain 9 / 16 residues).
- Files: results/phase8/showcase_multidomain_table.{csv,json,docx,pdf}; scripts showcase_multidomain_table.py, showcase_multidomain_pdf.js; s14_literature.py now also accepts a comma-separated ID list.

### Combined 21-protein table (2026-11-18, user request)
- Rows: last week's 12 (5 repository + 7 from the 60) + multi-domain showcase (6BYV, 1BBO, 2EW9, 5MS9, 1CRR, 2KUO, 2JT8, 6G81, 2V93); 2BW2 and 1B8T counted once.
- Reference mode run on the 8 new-scan entries without it (ref_mode_extra.py -> results/phase8/reference_mode_extra.csv; runs in results/expanded/runs_reference). 2V93: reference mode outputs two overlapping ranges (5..332, 111..370), reported as given. 6BYV reference result from the 60-protein study.
- 2KUO and 2JT8 labelled "RESUPER (partly)" vs real CYRANGE (coverage-only difference; unsupported extra split with a 1.87 Å range).
- Tally vs real CYRANGE: RESUPER 10, partly 4, CYRANGE 2, same 3, similar 1, unclear 1.
- Files: results/phase8/combined21_table.{csv,json,docx,pdf}; scripts combined21_table.py, combined21_pdf.js.

## S15: new repository (repo_v2, 2026-11-20) and publishable ~20-case table
- repo_v2 differs from the earlier repository only in src/resuper/cli.py and states.py (multiconformer
  state numbering/printing; tests). Domain detection unchanged. Check: v2 run with default parameters on all
  1,088 structures still on disk (222 random-sample + all differing scan entries): domains identical for
  1,088/1,088 (results/v2_check/v2_vs_v1.csv). Scan entries not on disk were identical to CYRANGE under v1;
  since the domain code is unchanged they are not re-run.
- Bug found and fixed in OUR wwPDB reader (scripts/wwpdb_cyrange.py): it read only model "1"; entries whose
  models are numbered otherwise (1CRR: 21-40) gave 0 CYRANGE domains. Affected: 1CRR, 1EZC, 1EZD only (none
  of the 222; the 108/108 validation is unaffected). Reader now uses the first model present. 1CRR, shown
  earlier as "CYRANGE reports no domain", was an artefact and is withdrawn. Evidence step re-run
  (previous outputs kept as *_before_modelfix.csv).
- Table rules (written after the S14 scan, before building the table; post hoc relative to the scan):
  A "RESUPER separates domains that CYRANGE merges": S14-qualifying; RESUPER reports more domains of
    >= 20 residues than CYRANGE; every RESUPER domain internally rigid (<= 1.0 A); no helix split;
    publication does not contradict.
  B "CYRANGE reports a non-rigid region as a domain": S14-qualifying; every RESUPER domain matches a
    CYRANGE domain (Jaccard >= 0.9); the remaining CYRANGE domain(s) have RMSD > 1.0 A and >= 8 residues;
    protein >= 60 residues; the original publication explicitly describes that region as disordered,
    flexible or unstructured. Candidates checked in order of decreasing extra-domain RMSD; entries with no
    relevant statement, a non-primary paper or an ambiguous statement are listed as not included.
### S15 outcome (2026-11-20)
- After the reader fix, 164 entries qualify (1CRR removed; 1EZC/1EZD do not qualify). Multi-domain count
  unchanged: 17 entries, 7 qualify (2JT8 excluded from A: RESUPER range 1.87 A), 10 not.
- Category B literature: included 9FLG, 1HN6, 2MCD, 7YWR, 1K9C, 2K8P, 2PXG, 2LWR, 2KTV, 5NR5. Not included
  (no explicit statement, unrelated or non-primary paper, or ambiguous): 2N47, 2LXD, 2K33, 2ECW, 1AP8, 2JY9,
  2RSM, 1N6Z, 2LHJ, 2BL5 (paper calls the C-terminal helix structured), 6X6N, 1HKO, 7BWE, 2K3N, 2N8X, 2KXR,
  2JVL, 2JQV, 1WFF, 1B8T, 2E0G, 1GXE, 2K4V, 2KHL, 1UFF, 1TWO, 2LF6, 2KHC, 6YI3, 5H3M, 1WJP, 9Y13.
- Table: results/phase9/resuper_vs_cyrange_cases.{csv,json,docx,pdf} (16 entries; scripts publishable_table.{py,js}).
- 2026-11-20 (user): table columns "Ensemble / annotation evidence" and "Original publication" merged into one
  "Evidence" column with numbered citations (Europe PMC metadata, results/phase9/refs/refs.json) and a
  reference list. Previous version kept as results/phase9/v1_resuper_vs_cyrange_cases.*.
- 2026-11-20 verification (scripts/verify16.py, results/phase9/verify16/): all 16 entries re-downloaded
  (PDB files and wwPDB validation reports into new directories, no cache), RESUPER repo_v2 re-run with
  default parameters, CYRANGE domains and RMSDs recomputed. 16/16: files identical to those used; RESUPER
  domains, CYRANGE domains, RMSDs and domain counts identical to the table; wwPDB domain counts consistent.
- 2026-11-20 (user): Methods + Results document remade as manuscript_methods/RESUPER_methods_results_JKMRS_v2.{docx,pdf}
  (script combined2_jkmrs.js): 16-case table (S15), 7SA5 multiconformer example (repo_v2 reporting; calibration
  caveat), selection process rewritten; 1A11 example and 1CRR removed; per-state result not included (user).
  Pin1 states re-checked with repo_v2: 2 states for 7SA5/7SUQ/7SUR/8SG2 (silhouette 0.745/0.744/0.667/0.622).
  Previous version kept (RESUPER_methods_results_JKMRS.{docx,pdf}).
- 2026-11-20 (user): v3 of the document: neutral wording of the illustrative-cases paragraph, Table 4 caption and Methods selection sentence; content and numbers unchanged. v2 kept.
- 2026-11-21 (user): repository package built: RESUPER_benchmark_cyrange_comparison.zip (benchmarks/cyrange_comparison/; paths made repository-relative; frozen-selection guard; Rmsd extra_dirs; verify16 uses shipped SHA-256 list; literature_check.csv instead of abstracts). Tested in a mock checkout: verify16 16/16 on all checks.

## S16: helix splits re-examined (2026-11-22, user request)
- Correction to S12/S14: a split inside one deposited HELIX record is not by itself an error. A
  continuous helix (same secondary structure) need not be one rigid body; the user checked 2KDR and its
  two segments do move relative to each other. "Single helix" statements in the 2KDR and 9G7T papers do
  not address the relative motion of the segments, so they do not contradict RESUPER.
- Check (scripts/helix_split_motion.py -> results/expanded/helix_split_motion.csv): all 21 helix-split
  entries of the 908, 24 junctions. Relative-orientation scatter of the second segment after fitting the
  first exceeds that of the two halves of the longer segment at every junction (ratio 1.26-44; 2KDR 2.3,
  9G7T 2.1). Helix H-bonds across the junction stay intact (>= 0.9 of conformers) at 12 junctions.
- The Table 4 rule (no split inside one HELIX record) is kept as a conservative selection rule; the
  manuscript Discussion no longer calls these splits over-splitting.
- results/phase9/literature_check.csv: 2KDR and 9G7T relabelled "not contradicted (S16)".
