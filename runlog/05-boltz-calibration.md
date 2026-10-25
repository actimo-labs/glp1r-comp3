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

---

# ADDENDUM — the shuffle null gave the wrong answer

The folded-non-binder control (6 comp-02 TNF designs scored against GLP1R,
`controls/boltz_negbind/`) **reverses the conclusion above.**

| Metric | POS native | classC best | folded non-binder | discriminates? |
|---|---|---|---|---|
| `pae_interaction` | 4.56 | 12.55 | **20.04 – 28.66** | **yes** |
| `ipTM` | 0.949 | 0.942 | **0.113 – 0.561** | **yes** |
| `plddt_binder` | 91.6 | 84.1 | **75.9 – 86.8** | **no, overlaps** |

`plddt_binder` was measuring foldedness. Our class C designs at 84.1 sit inside
the folded-non-binder band (median 81.5) and cannot be separated from it. The
shuffle made `plddt_binder` look like the best metric purely because a shuffled
sequence does not fold — it was answering the wrong question.

`pae_interaction`, which the shuffle made look useless, is the real
discriminator: every folded non-binder scores ≥ 20.04, our class C designs
12.55–13.85.

**Final cuts**: primary `pae_interaction < 20.0`; secondary `ipTM > 0.56`;
`plddt_binder > 70` retained only as a did-it-fold sanity check.

**On using ipTM.** Brief §3 says rank on PAE and pLDDT "not ipTM", and comp-02's
non-negotiables forbid global ipTM outright. Both were reasoning about a
homotrimer, where global ipTM is dominated by target-internal chain pairs
sitting at ~0.9 and is therefore meaningless. GLP1R is a single chain: there is
no target-internal pair, so ipTM **is** the binder–target score. The separation
is 0.94 vs ≤0.56. Used as a secondary, with this reasoning stated rather than
the rule quietly ignored.

**Class B re-read.** `classB_0_s1` scores pae_interaction 19.36 and ipTM 0.653 —
passing both cuts, marginally — where AF2-ig gave it a flat, uninformative 27.4.
Further evidence that AF2-ig should not gate the TMD classes.
