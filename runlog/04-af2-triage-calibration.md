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
