# RESUPER vs CYRANGE: benchmark on PDB NMR ensembles

This folder holds the scripts, frozen inputs, decision rules and results behind the comparison of
RESUPER (default method) with CYRANGE in the manuscript: the random-sample tests (Table 3) and the
illustrative cases (Table 4). It only **reads** the RESUPER code in `../../src`; nothing here changes it.

## What was compared

| Comparator | Source |
|---|---|
| RESUPER (default) | this repository, `python -m resuper <file>.pdb`, default parameters |
| CYRANGE | per-residue domains (`cyrange_domain_id`) in each entry's **wwPDB NMR validation report**; read from the first model present. Check: identical to all 108 CYRANGE 2.0 reports in `testdata/` |
| Reference mode | `--method reference`, an approximate reimplementation of CYRANGE 2.0 (comparison only; its results are not CYRANGE's) |

## Main results (already in `selection_study/results/`)

- **Random samples (222 structures; primary test set 162)**: RESUPER and CYRANGE identical for 136/222;
  no significant difference on the four prespecified hypotheses (specified with the reference mode as comparator; post hoc for CYRANGE) (annotated domains recovered,
  confirmation sample: RESUPER 32, CYRANGE 38, p = 0.125). `real_cyrange/tests.txt`
- **Illustrative cases (16)**: from 8,709 further entries, cases where independent evidence supports
  RESUPER, in two categories (A: RESUPER separates domains CYRANGE merges; B: CYRANGE reports a region
  described as disordered as an extra domain). `phase9/resuper_vs_cyrange_cases.{csv,pdf,docx}`
- **Version check**: the current repository gives identical domains to the earlier version on all
  1,088 structures re-run. `version_check/v2_vs_v1.csv`; all 16 cases re-verified from fresh downloads,
  `phase9/verify16/verify16.csv`

These cases illustrate where RESUPER performs better; the random samples estimate how often.

## Layout

```
unpack_runs.py                unpacks the bundled run outputs (*.zip) in place
scripts/                      shared helpers (paths, report parser, comparison, structure tools)
data/                         manifest of the earlier 37-protein benchmark (motion baseline)
results/phase4_runs, phase5_comparison   earlier benchmark outputs used for the motion baseline
selection_study/
  docs/decisions_log.md       every rule, written down before it was applied (S1-S15), and deviations
  docs/real_cyrange_results.md
  scripts/                    pipeline (below)
  results/                    frozen lists and all results used in the paper
```

Downloaded PDB files, validation reports and RCSB annotations are not included (they are re-downloaded
by the scripts). Publication abstracts are not redistributed; `phase9/literature_check.csv` lists, for
every candidate checked, the PMID, the decision and the statement or reason.

## Re-running

Requirements: Python ≥ 3.10, numpy, scipy (`pip install -r requirements.txt`), internet access to
files.rcsb.org, data.rcsb.org and www.ebi.ac.uk. Run everything from this folder:

```bash
cd benchmarks/cyrange_comparison
```

Quick check (about 5 min): re-download the 16 table entries, re-run RESUPER and CYRANGE, compare with the table:

```bash
python selection_study/scripts/verify16.py
```

Full pipeline: first unpack the bundled run outputs (kept as zip archives so this folder stays under
GitHub's 100-files-per-upload limit on the website):

```bash
python unpack_runs.py
```

Then, in order (frozen lists are used where they exist, so selections are not redrawn):

| Step | Command | Output |
|---|---|---|
| Candidate pool and first sample (seed 20261111) | `python selection_study/scripts/phase3_select_candidates.py` | `results/phase3/` |
| Annotations | `python selection_study/scripts/phase3_annotate.py` | `results/phase3/` |
| Screening and runs, first sample | `python selection_study/scripts/phase4_screen_and_run.py` | `results/phase4/` |
| Comparison, first sample | `python selection_study/scripts/phase5_compare.py` | `results/phase5/` |
| Confirmation sample (seed 20261113) | `python selection_study/scripts/confirm_select.py` then `confirm_run.py` | `results/confirm/` |
| CYRANGE from wwPDB, validation vs 108 saved reports | `python selection_study/scripts/wwpdb_cyrange.py validate` | `results/wwpdb_cyrange/` |
| Random-sample tests vs CYRANGE (Table 3) | `python selection_study/scripts/real_cyrange_all.py` | `results/real_cyrange/` |
| Scan of the remaining pool (~2 h, resumable) | `python selection_study/scripts/expanded_scan.py --workers 6` | `results/expanded/scan.csv` |
| Evidence and ranking (~25 min) | `python selection_study/scripts/expanded_evidence.py` | `results/expanded/{evidence,ranked}.csv` |
| Literature for candidates | `python selection_study/scripts/s14_literature.py ID1,ID2,...` | `results/expanded/literature/` |
| Table 4 data / Word table | `python selection_study/scripts/publishable_table.py`; `NODE_PATH=$(npm root -g) node selection_study/scripts/publishable_table.js` | `results/phase9/` |
| Version check (needs the structures downloaded by the steps above) | `python selection_study/scripts/version_check.py` | `results/version_check/` |

`RESUPER_REPO=/path/to/RESUPER` overrides the repository location if this folder is moved.

## Notes and limitations

- The four hypotheses were registered with the reference mode as comparator; applying them with
  CYRANGE (wwPDB) uses identical rules but is post hoc (`decisions_log.md`, S11, S13).
- Category B of Table 4 was defined after the scan and before the table was built (S15); inclusion
  required an explicit statement in the original publication.
- A reader error (only model "1" read) affected 3 entries (1CRR, 1EZC, 1EZD) and was fixed before the
  final table (S15). It did not affect the random samples or the 108-report validation.
- Domain RMSDs for both methods are computed with the same routine; a lower RMSD alone does not mean a
  better assignment.
- `decisions_log.md` is the original working log; paths in it refer to the working environment.
