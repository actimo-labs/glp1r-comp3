# Run log 06 — Submission spec, confirmed from the competition page

Fetched 2026-10-25 from
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/glp1r
This closes one of the three items brief §1 left open.

## Required CSV columns (FAQ 5)

Ordered by rank, **top row = highest ranked**. Minimum columns:

| Column | Value for us |
|---|---|
| `name` | unique identifier |
| `sequence` | the designed protein sequence |
| `molecule_class` | `single_chain` (options: single_chain, nanobody, scfv, fab_kappa, fab_lambda) |

`molecule_class` is an optional column and can also be set per design after
upload. scFvs and Fabs are submitted as a single sequence written `{VH}:{VL}`;
ours are all single-chain proteins so no colon appears in any sequence.

**Correction, 2026-10-26.** An earlier version of this file recorded the first
accepted value as `protein`. That was wrong — it came from a paraphrased fetch
of the FAQ rather than the upload form itself, and `protein` is not an accepted
value. The correct value is `single_chain`, confirmed by Harjas against the
form. All 20 designs were rewritten accordingly.

Methods, metrics and structure models are optional but encouraged — we supply them.

## Constraints

- length **10–250 aa**, single chain
- **up to 20** designs for Track 3
- de novo only; see the FAQ 8 ruling in `runlog/07-denovo-ruling.md`
- binding to a tag does not count; the FLAG and His tags on the receptor
  construct, and the design's own C-terminal Twin-Strep, are off limits as
  epitopes

## Twin-Strep tag — DECIDED: submit WITHOUT it

The page refers to "your design's C-terminal Twin-Strep capture tag" but never
states whether entrants include it in the submitted sequence. Harjas decided
2026-10-25: **submit the designed sequence only, no tag.**

Reading behind it: the tag is assay infrastructure for immobilisation, appended
by the organisers to every entry uniformly. Including it would also consume
length budget against the 250 aa cap for no benefit.

Consequence that is already handled: every design must tolerate a C-terminal
fusion, so nothing may depend on a free C-terminus (brief §4: "Any design whose
mechanism needs a free C-terminus is dead on arrival"). Our designs place the
functional element at the **N**-terminus, so this holds.

## Assays (for the methods report)

- **Binding**: design immobilised, receptor in solution, pH 7.4. Counts as
  binding if the sensorgram gives a fittable KD, or the association shift
  clears 300% over negative control.
- **Activation**: cell-based CRE-SEAP reporter, cAMP-driven, on confirmed
  binders only, versus untreated and a reference agonist.

## Ranking

1. lowest KD among all binders, agonist or not; groove binding counts in favour
2. activation relative to the reference agonist
3. winner must do both: "the stronger activation wins, then the higher affinity"
