# Run log 07 — De novo ruling on class B, and the rebuild

## The rule

Competition FAQ 8, fetched 2026-10-25:

> "They must be de novo and zero-shot, defined as **not having relied on an
> initial binder as a starting point** ... having adequate **sequence- and
> structural-diversity** [from known proteins] ... You may not take an existing
> binder and modify it. Designs must be produced from scratch."

and, permitting the opposite use:

> "you may take known binders and use them to **evaluate/calibrate design
> filtering metrics**" and "to fine-tune/train models."

## What we had built

Class B v1 held **GLP-1 residues 7–15 backbone coordinates fixed** as an
RFdiffusion motif, with all sequence designed de novo. Brief §5 had decided to
"design the message module de novo at the pharmacophore level — reproduce the
*interactions*, not the sequence", and a backbone graft goes beyond that.

## Ruling: non-compliant, withdrawn

The deciding clause is "sequence- **and structural**-diversity". Class B v1 had
full sequence diversity and deliberately **reduced structural diversity at the
functional site** — 9 residues of backbone taken from GLP-1, which is the native
binder of this exact target. The permissive reading requires arguing GLP-1 is
"not a binder", which does not survive contact with the rule.

Note the asymmetry FAQ 8 draws, which we are on the right side of elsewhere:
known binders **are** sanctioned for calibrating filter metrics. That is exactly
how `controls/negbinders.csv` uses comp-02's TNF designs, and how the native
GLP-1 positive control is used. Those stay.

## The rebuild — and it is better

`scripts/32_generate_classB.sh` uses **no ligand coordinates at all**.
RFdiffusion is pointed at the deep activation-trigger residues on the receptor
itself and must find its own way in:

```
A387 (35.1 A)  the free-alpha-amine partner, the activation trigger
A234 (29.3 A)  direct peptide contact
A233 (29.0 A)  deep pocket
A391 (38.4 A)  deep
A241 (39.5 A)  deepest
```

Everything derives from 6X18 chain R, which the competition lists as the
structure reference.

### Pilot, 3 designs

| Design | len | contacts | deepest | triggers hit | lipid frac | QC |
|---|---|---|---|---|---|---|
| classBL90_0 | 90 | 13 | **43.7 Å** | **3** (241, 387, 391) | 0.08 | ok |
| classBL90_1 | 90 | 10 | **44.4 Å** | **3** (241, 387, 391) | 0.00 | ok |
| classBL130_0 | 130 | 10 | 40.6 Å | 1 (241) | 0.20 | ok |

Against the withdrawn motif version: deepest **30.4 Å**, **0** trigger residues.

E387 is contacted at **4.12 Å** and **4.40 Å**. The compliant formulation
reaches the activation trigger and the motif version did not — the fixed motif
appears to have been *constraining* the search rather than guiding it.

Length 130 is worse on both trigger count and lipid contact. Production uses
**90 and 110**.

## Status

Class B keeps its 10 slots. Production queued behind C+A so it runs
uncontended.
