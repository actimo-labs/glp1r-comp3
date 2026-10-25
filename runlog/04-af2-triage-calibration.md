# Run log 04 — AF2-ig triage, calibrated on controls

48 designs (32 class C, 6 class B, 2 positive, 8 negative controls), 8 distinct
lengths, 1228 s total. The length-sorting fix worked: 8 compiles rather than 48.

## The missing FastRelax has a measured cost

Native GLP-1, same sequence, same backbone, the only difference being whether
binder side chains are present:

| Representation | pae_interaction | plddt_binder |
|---|---|---|
| full atom (native structure) | **7.06** | 84.0 |
| backbone only (what our designs are) | **17.39** | 51.5 |

Upstream dl_binder_design runs MPNN → **Rosetta FastRelax** → AF2-ig, and
FastRelax builds the side chains. We have no PyRosetta, so every design reaches
AF2 as backbone + sequence only. This is the price, and it is large.

The practical consequence: **no published threshold transfers.** Brief §3 quotes
Baker's pae_interaction < 8–10. In our representation the native agonist itself
scores 17.39, so that cutoff rejects a validated full agonist. Thresholds come
from the controls or they are meaningless.

## Calibrated thresholds

| Control | pae_interaction |
|---|---|
| POS native GLP-1 | **17.39** |
| POS exendin-4 (on GLP-1's backbone) | 27.33 — see caveat |
| NEG shuffles (n=8) | min **22.00**, median 23.78 |

The exendin-4 positive failed. That is expected rather than alarming: exendin-4
was threaded onto GLP-1's backbone, which is not its own, so it is a
non-native pairing and a weak control. Treat the GLP-1 positive as the only
valid ceiling. One usable positive is thin, and the methods report should say so.

**Cut: pae_interaction < 22.00** — below every negative.

## Results by class

| Class | n | best pae_int | median | verdict |
|---|---|---|---|---|
| C (ECD groove) | 32 | **17.52** | 26.69 | **real signal** |
| B (two-domain) | 6 | 27.42 | 27.58 | no signal |
| A (vestibule rim) | 20 | 27.26 | 27.9 | no signal |

**7/32 class C designs beat every negative**, best at 17.52 against the positive
control's 17.39 — indistinguishable from the native agonist in the same
representation. classC_3 contributed 4 of the 7.

## Why A and B show nothing — and why it is not evidence the designs are bad

AF2's PAE saturates near 31.75. Scores of 27–28 are not "weak interface", they
are "no information". Every class A and class B design sits there.

The split is explained by *where* the epitope is, not by design quality:

* class C targets the **ECD**, a soluble globular domain — AF2's home ground.
* classes A and B target the **TMD** — and brief §3's own benchmark puts AF2-ig
  at median DockQ **0.03** on 113 GPCR complexes, against Boltz-2 at 0.56.

So this is AF2's documented blind spot on membrane-protein interfaces. Using it
to reject class A and B designs would be reading a known-broken instrument.

## Consequence for the pipeline

- **class C**: AF2-ig is a valid gate. `pae_interaction < 22.0`.
- **classes A and B**: AF2-ig is **not** a gate. Keep `plddt_binder` only, as a
  check that the binder folds at all, and let Boltz-2 score the interface.
- That moves A and B's full load onto Boltz-2, which is more expensive. The
  pre-filter for those classes becomes interface QC + MPNN score + plddt_binder,
  then Boltz on the top slice.

Recorded in the "not computed / not run" list: no FastRelax; its measured cost
on the positive control is pae_interaction 7.06 → 17.39.

---

# ADDENDUM — AF2 checked against the folded-non-binder null too

Three comp-02 TNF designs threaded onto class C backbones of matching length
and run through AF2-ig:

| Control | pae_interaction | plddt_binder |
|---|---|---|
| NEGBIND_tnf_186_s4 | 26.10 | 92.3 |
| NEGBIND_tnf_6_s3 | 26.86 | 94.1 |
| NEGBIND_tnf_214_s3 | 27.45 | 94.3 |

Two things follow.

**The `pae_interaction < 22.0` cut survives.** Folded non-binders score 26.1–27.45,
comfortably above it, and above the shuffle null's 22.00 as well. The cut sits
below both nulls, so it is conservative in the right direction. Kept unchanged.
It could be relaxed to ~24 to admit more class C candidates, but all seven class
C designs that currently pass score 17.52–21.90, so relaxing buys nothing and
costs Boltz time.

**It confirms AF2-ig is blind on the TMD classes, now against a proper control.**

```
folded non-binders (definitely do not bind) : 26.10 - 27.45
our class A and class B designs  (n = 26)   : 27.26 - 28.35
```

Our TMD designs are **indistinguishable from proteins known not to bind** —
in fact marginally worse. Earlier this was inferred from PAE saturation; it is
now measured against a real negative. Gating class A or B on AF2-ig would be
rejecting designs on a signal the instrument demonstrably cannot produce for
this epitope. Boltz-2 does that job.

Note also `plddt_binder` 92–94 for the non-binders, higher than any of our
designs. Same lesson as the Boltz calibration: `plddt_binder` tracks foldedness,
not binding, in both models.

---

# ADDENDUM 2 — scoring class C against the full receptor produced false negatives

Running production triage with the binder threaded onto the **full 384-residue
receptor**, all 197 scored class C designs failed the cut, best 26.60. The
pilot, which threaded onto the **isolated 101-residue ECD**, had 7 of 32 pass
with a best of 17.52. Something was wrong with the representation, not the
designs.

Diagnostic: the same 8 designs, both ways.

| Design | full receptor | ECD only |
|---|---|---|
| classCL65_35_s2 | 27.03 | **6.39** |
| classCL65_4_s1 | 27.08 | **8.58** |
| classCL65_21_s3 | 27.09 | 25.46 |
| classCL65_4_s2 | 26.60 | 26.42 |
| classCL65_26_s3 | 26.76 | 26.71 |
| classCL65_65_s2 | 27.05 | 26.80 |
| classCL65_4_s4 | 26.88 | 26.90 |
| classCL65_33_s3 | 27.08 | 27.48 |

Two designs scoring **better than the native GLP-1 positive control (17.39)**
were being reported as complete failures. Adding 283 residues of transmembrane
domain — which is not part of the class C epitope — pushes the whole complex
into AF2's GPCR blind spot and the interface signal disappears.

Note how the full-receptor numbers cluster in 26.60–27.09. That is a noise band,
not a ranking. The 8 designs above were chosen as the "least bad" by
full-receptor score, i.e. effectively at random, and 2 of 8 turned out to be
excellent — which implies many more good designs were buried in the other 604.

**Fix.** Class C AF2 triage is run against the **ECD alone**, which is both the
representation AF2 can handle and the honest target: the class C epitope lies
entirely within the ECD. Boltz-2 still validates against the complete receptor,
so TMD clashes would still be caught.

**Generalisation.** The target representation used for scoring must match the
epitope, not default to "the biggest structure available". Classes A and B
genuinely need the TMD and are therefore not AF2-scorable at all; class C does
not, and is.
