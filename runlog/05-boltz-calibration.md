# Run log 05 — Boltz-2 calibration

5-design pilot, receptor pinned to the active-state template
(`target/glp1r_R_full.cif`, 6X18 chain R). **130.4 s total, 26.1 s/design.**

That is ~4x faster than the 100 s/design I had budgeted, which is why classes A
and B could be moved onto Boltz-2 without breaking the schedule.

First attempt failed instantly: `FileNotFoundError: 'boltz'`. The CLI is not on
PATH inside a subprocess even when the script is run by that env's python. Now
invoked by absolute path.

## Results

| Design | pae_interaction | plddt_binder | ipTM |
|---|---|---|---|
| **POS** native GLP-1 | 4.56 | **91.6** | 0.949 |
| **NEG** shuffled GLP-1 | 6.76 | **47.1** | 0.898 |
| classC_3_s1 | 12.55 | 84.1 | 0.942 |
| classC_3_s4 | 13.85 | 84.8 | 0.939 |
| classB_0_s1 | 19.36 | 61.8 | 0.653 |

## pLDDT_binder discriminates; PAE does not

Brief §3 Stage 4 says to score on "inter-chain PAE and peptide/binder pLDDT".
Measured here, those two are not equally useful:

* **pae_interaction**: positive 4.56 vs negative 6.76 — a 2.2 Å gap, and the
  class C designs (12.5–13.8) sit *further* from the positive than the negative
  does. Ranking on PAE would put a shuffled peptide above every real design.
* **plddt_binder**: positive 91.6 vs negative 47.1 — a 44-point gap, with the
  class C designs at ~84, clearly on the positive side.
* **ipTM**: 0.949 vs 0.898. Useless, exactly as the brief predicted.

This is the brief's own caveat showing up in our data: Boltz-2 confidence
"systematically overestimate[s] correctness for misplaced peptides", and PAE in
particular cannot tell medium from high quality. **Ranking weight goes on
plddt_binder.** PAE is recorded and used only as a weak secondary.

## Known weakness in the null, and the fix

A shuffled sequence does not fold, so `plddt_binder` on it is measuring
foldedness, not binding. A well-folded non-binder would score high and pass.
The shuffle calibrates the wrong failure mode.

Better null: a real de novo designed miniprotein with no reason to bind GLP1R.
comp-02 has 1640 such sequences on disk — TNF binders, folded, designed by the
same ProteinMPNN install. Six sampled (60–90 aa) and scored against GLP1R as
`controls/negbinders.csv`. Whatever `plddt_binder` those earn is the floor for
"folds well, does not bind this target", which is the threshold that actually
matters.

## Thresholds (to be finalised once the folded-non-binder null returns)

- primary: `plddt_binder`, cut above the folded-non-binder distribution
- secondary: `pae_interaction`, recorded, weak
- `ipTM`: recorded, never used to rank — consistent with the comp-02
  non-negotiable about not ranking on global ipTM
