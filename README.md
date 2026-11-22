# RESUPER

RESUPER identifies the ordered regions (rigid domains) of a protein structure
ensemble, selects the residue range to superimpose for each domain, and detects
distinct conformational states. It reads multi-model PDB or mmCIF files: NMR
bundles, molecular-dynamics snapshots or ensemble-refinement models.

```
    Optimal range 5..39: RMSD 0.92 A, 0 gaps, 35 residues.
    Optimal range 53..163: RMSD 0.87 A, 0 gaps, 111 residues.

    Conformational states: 2 (populations 0.50, 0.50); conformers 2, 4-5, 8-9, 12-13, 16, 18, 20 -> state 1, conformers 1, 3, 6-7, 10-11, 14-15, 17, 19 -> state 2.
```

Each range is a set of residues that move together and can be fitted by one
rigid-body superposition. Flexible termini, linkers and hinges are left out, so
the reported RMSD describes the well-defined part of the structure.

## Graphical interface in POKY

RESUPER is distributed with POKY, which provides its graphical interface. Download
POKY from https://poky.clas.ucdenver.edu and open the RESUPER dialog (two-letter code yS). 
From the dialog you can set the parameters, detect conformational states, run the analysis
per state, save the superposed ensemble (one file per state if you wish), show
the ensemble colored by state in PokyMol, and copy the ranges formatted for
Xplor-NIH, CYANA and PyMOL.

## Installing the command-line program

Python 3.9 or later with NumPy and SciPy is required.

```bash
git clone https://github.com/pokynmr/RESUPER_deploy.git
cd RESUPER_deploy
pip install .
resuper --help
```

To run from the source tree without installing:

```bash
PYTHONPATH=src python3 -m resuper ensemble.pdb
```

To run the test suite:

```bash
pip install pytest
PYTHONPATH=src python3 -m pytest tests -q
```

## Usage

```
resuper [PARAMETERS] <ensemble.pdb | ensemble.cif>
```

```bash
resuper ensemble.pdb                                 # ordered regions and ranges
resuper --states ensemble.pdb                        # also count conformational states
resuper --per-state ensemble.pdb                     # ranges of each state (implies --states)
resuper --write-ensemble superposed.pdb ensemble.pdb # write the superposed ensemble
resuper --split-states states/ --states ensemble.pdb # one superposed file per state
```

The parameters use CYRANGE's names and units:

| option | meaning | default |
|---|---|---|
| `-a` | absolute RMSD decrease (δ_abs) | 1.6 Å |
| `-b` | domain-boundary extension (m) | 3 residues |
| `-d` | relative-decrease base (δ) | 1.2 |
| `-g` | gap-opening penalty (γ) | 0.4 |
| `-u` | RMSD threshold that stops refinement | 0.1 Å |
| `-gw` | gaps shorter than this are refilled (g) | 3 residues |
| `-m` | minimum domain size (μ) | 8 |
| `-r` | restrict the analysis to a residue range | all |
| `-l` | custom torsion library | internal |

`--method reference` is an approximate reimplementation of CYRANGE 2.0, kept for
side-by-side comparison only. Its results are not those of CYRANGE.

## Method in brief

1. Core atoms: Cα atoms of residues with a torsion whose angular order parameter
   reaches the Q_i-optimal cutoff of Kirchner and Güntert (BMC Bioinformatics 12,
   170, 2011).
2. Domains: the core atoms are clustered by Ward linkage on the variances of their
   distances. A node of the tree is split when the split score
   r(node) − r(left) − r(right) is at least 0.12 Å, when the distance-variance
   matrix shows two blocks in relative motion, or when an acceptable split exists
   further down the tree.
3. Ranges: each domain starts from its core span extended by three residues and
   is refined by the published gap-weighted RMSD-decrease loop on N, Cα and C′.
4. States: the conformers are clustered on their Cα–Cα distance patterns; the
   number of states is the partition with the largest silhouette increase,
   reported when its silhouette is at least 0.60.

## Repository layout

- `src/resuper/` – the package (`core.py` algorithm, `cli.py` command line,
  `states.py` state detection, `writer.py` superposed output).
- `tests/` – test suite.
- `testdata/` – ensembles and captured reports used by the tests and benchmarks,
  including the Pin1 ensembles 7SA5, 7SUQ, 7SUR and 8SG2.
- `scripts/` – benchmarks with known domains and the state analyses.
- `benchmarks/cyrange_comparison/` – comparison with the CYRANGE results in the
  wwPDB validation reports (random samples and illustrative cases).

## Benchmarks

Run every command from the repository root. Some inputs are not stored in this
repository and are downloaded first:

```bash
python3 scripts/pacsy_eval.py --fetch            # PACSY-selected PDB entries (testdata/pacsy_ids.txt), from the RCSB PDB
python3 scripts/aisar_bundles.py --fetch --testdata   # AISAR two-state bundles
```

| result | commands |
|---|---|
| 134 real NMR ensembles | `python3 scripts/eval_methods.py` |
| held-out sets 1 and 2, test set 3 | `python3 scripts/heldout/heldout.py build [--set 2\|3]`, then `... score [--set 2\|3]` |
| test set R (motions applied to real ensembles) | `python3 scripts/heldout/realhinge.py build`, then `python3 scripts/heldout/realhinge_score.py` |
| multi-body set (788 ensembles) | `python3 scripts/bench/bench.py build-base`, `... build`, `... run`, then `... score` |
| perturbation stability | `python3 scripts/heldout/perturb.py build`, then `... score` |
| sensitivity to τ, linkage and torsion set | `python3 scripts/stagea/sens_split.py` |
| 13 synthetic ensembles with known domains | `python3 scripts/synthetic/truthbench.py --capture`, then `python3 scripts/synthetic/truthbench.py` |
| Pin1 domain separations and per-state RMSD | `python3 scripts/states_motion.py testdata/7sa5.pdb`, `python3 scripts/states_rmsd.py testdata/7sa5.pdb` |
| random samples and illustrative cases | see `benchmarks/cyrange_comparison/README.md` |

Generated ensembles are written to `cywork/`. The CYRANGE reports that the
benchmarks compare against are included (`testdata/*.out`, `testdata/oracle_*`),
so scoring does not need CYRANGE. Commands named `capture` (and
`truthbench.py --capture`, which builds the synthetic ensembles and captures their
reports in one step) run a licensed CYRANGE binary on the host set in
`RESUPER_CYRANGE_HOST`.

The test suite skips the tests whose inputs have not been downloaded.

## Citation

If you use POKY's RESUPER GUI, please cite POKY:

W. Lee, M. Rahimi, Y. Lee, and A. Chiu, Bioinformatics 37, 3041 (2021)

A paper describing RESUPER is in preparation.

## License

RESUPER is part of POKY and is available under the POKY dual license:

- academic and non-commercial use: POKY Academic License (`LICENSE.POKY.ACADEMIC`);
- commercial use: POKY Commercial License (`LICENSE.POKY.COMMERCIAL`). Contact
  cuinnovations@ucdenver.edu and reference POKY software or CU5598D.

Copyright © 2021-present, The Regents of the University of Colorado.
