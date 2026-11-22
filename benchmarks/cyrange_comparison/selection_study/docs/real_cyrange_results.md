# RESUPER vs real CYRANGE 2.0 on all 222 proteins (decision S13)

Script: `scripts/real_cyrange_all.py`. Outputs: `results/real_cyrange/`
(per_protein.csv, method_specific_domains.csv, wwpdb_status.csv, tests.txt).

## Source of real CYRANGE results
The wwPDB validation report of every NMR entry stores CYRANGE's domain for each
residue. Check (S12): for all 108 PDB entries with a CYRANGE 2.0 report saved in
the RESUPER repository, the wwPDB domains are exactly identical (9/9 multi-domain).
All 222 proteins (60 original + 162 confirmation) have a result; no residue was
lost when mapping to our files. Domain RMSDs were computed here (reference-mode
atom convention), not taken from CYRANGE.

## Status of the tests
H1-H4 were pre-registered against reference mode (S11). Running them against
real CYRANGE uses the same rules but is a post hoc analysis. The 162-protein
sample is the cleaner test (the 60 were looked at in depth earlier).

## Results: confirmation sample (162)
| | RESUPER vs CYRANGE |
|---|---|
| Identical / minor / substantial difference | 99 / 33 / 30 |
| H1 annotated domains recovered (A+B, 47 proteins) | 32 vs 38; RESUPER better 0, CYRANGE better 4; p = 0.125 |
| H2 annotated domains missed | 6 vs 4; p = 0.5 |
| H3 over-splitting in controls (C+E, 87) | RESUPER only 3, CYRANGE only 4; p = 1 |
| H4 short helix reported (D, 28) | RESUPER only 0, CYRANGE only 2; p = 0.5 |
| Domains reported by only one method | RESUPER 8 (5 pass two-part check); CYRANGE 16 (11 pass) |

Reference mode vs real CYRANGE on the same 162: identical 117, substantial 34.

## Results: all 222
- RESUPER identical to real CYRANGE for 136 (61%), same domain count for 196 (88%).
- Reference mode identical to real CYRANGE for 158 (71%): reference mode is
  closer to CYRANGE than RESUPER is, but differs substantially for 52.
- Annotated domains recovered (A+B, 71): RESUPER 63, CYRANGE 66; RESUPER better 4,
  CYRANGE better 6; p = 0.68. None of H1-H4 is significant.
- Domains reported only by one method: RESUPER 16 (7 pass the two-part check),
  CYRANGE 29 (18 pass).
- Group E (peptides): RESUPER identical to CYRANGE for only 16 of 55.

## What this changes
1. No general advantage for RESUPER over real CYRANGE; if anything a small,
   non-significant tendency the other way (more annotated domains recovered by
   CYRANGE in the confirmation sample: 4 proteins vs 0).
2. The earlier descriptive pattern (RESUPER's method-specific domains look like
   genuine rigid units more often) reverses: against real CYRANGE, CYRANGE's own
   method-specific domains pass the two-part check more often (18/29 vs 7/16).
3. Several earlier results came from reference mode's imitation errors, e.g. its
   over-splitting of controls (group C: reference mode 64 domains, CYRANGE 57,
   RESUPER 57) and its merging of multi-domain proteins (group B).
4. Proteins where RESUPER recovers more annotated domains than real CYRANGE: 2BW2,
   5JS8, 1WHE (RESUPER input problem; not assessable), 6BYV. The reverse: 2L9N,
   2RPC, 2DKT, 2DLQ, 1Z87, 2KUE.
