# Protocol: making the RESUPER default method state of the art

Written before the improvement loop started (2026-10-21). The comparison is with CYRANGE,
the established tool for this task (also used by the wwPDB NMR validation pipeline); no
other residue-range program is available to run.

## Data roles

- **Development** (may be looked at and tuned on): the ubiquitin truth bench
  (`scripts/synthetic/truthbench.py`), held-out sets 1 and 2 (`scripts/heldout/heldout.py`),
  the 134 real ensembles (`scripts/eval_methods.py`) and their perturbations
  (`scripts/heldout/perturb.py`). Sets 1 and 2 stopped being held out once a method change
  was derived from them.
- **Test**: a new held-out set is generated, and its design committed, after the
  candidate method is frozen and before the candidate is run on it. A test set is used
  for one verdict only; if the verdict fails and the method is revised, the old test set
  becomes development data and a new test set is made.

## State-of-the-art criteria (all must hold on the test set)

On ground-truth ensembles (96 scenario ensembles plus the 56-ensemble hinge sweep):
1. domain count correct: RESUPER >= CYRANGE;
2. mixed domains: RESUPER <= CYRANGE;
3. contamination (flexible residues selected): RESUPER <= CYRANGE + 0.005;
4. mean Jaccard with the true domains: RESUPER >= CYRANGE;
5. hinge sweep: number of two-domain detections summed over all amplitudes: RESUPER >=
   CYRANGE.

On the 134 real ensembles (development data, checked for every candidate):
6. same domain count as CYRANGE on at least 130;
7. no catastrophe: no ensemble whose residue-weighted domain RMSD exceeds CYRANGE's by
   more than 1.0 A;
8. under the five perturbations, the number of ensembles whose domain count changes is
   at most CYRANGE's.

And: the reference mode is unchanged and pytest passes.

## Round 2: benchmark B ("beat CYRANGE all around", 2026-10-21)

Development data now also includes benchmark B (`scripts/bench/bench.py`, 40 PDB
proteins, v2 with side-chain disorder). The test is benchmark C: the same design on 40
further PDB proteins drawn with a new seed and excluding benchmark B's pool and
testdata/, generated and committed after the candidate is frozen and before it is run.

Criteria on benchmark C, for every suite (clusters, expclus, idr, idp, multi), against
the CYRANGE binary and against RESUPER's reference mode:
1. domain count correct >=;
2. mixed domains <=;
3. adjusted Rand index >=;
4. contamination <= (+0.005);
5. Jaccard >=.
Errors suite: no crash, and at least as many acceptable outcomes (near for defects with
a sensible answer; any non-crash outcome otherwise) as the CYRANGE binary.
Plus the round-1 criteria 6-8 on the real ensembles, the reference mode unchanged,
and pytest passing.
