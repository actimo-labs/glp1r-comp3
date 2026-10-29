The pipeline, in five steps

1. Prepare the target. Took the cryo-EM structure of GLP1R (PDB 6X18) with its natural hormone GLP-1 bound. Worked out where on the receptor to aim by measuring which receptor residues GLP-1 actually touches — computed from the structure rather than copied from papers.

2. Generate backbones. RFdiffusion invents protein shapes that fit against a chosen patch of the target. We made 568 backbones across three strategies.

3. Design sequences. ProteinMPNN works out which amino acids would fold into each shape. ~1,900 sequences.

4. Filter hard. Several independent checks, each calibrated against controls (below). Most designs die here.

5. Rank and cut to 20. Best survivors, with diversity rules so we don't submit 20 variations of one design.

---

The decisions that mattered

Where to aim. GLP1R has two parts: an outer domain (easy to bind, like a landing pad) and a deep pocket buried in the membrane (hard to reach, but it's where activation actually happens). We split the portfolio to cover both — 10 designs on each.

We validated against controls, not published thresholds. We ran the real GLP-eline as a "known good", and ran unrelated real proteins through as "knownbad". This mattered enormously: the standard published cutoff would have rejected the natural hormone itself. Any threshold we'd borrowed would have been wrong.

We abandoned a whole strategy mid-run. One approach (250 designs' worth) used GLP-1's own backbone as a template. Reading the competition rules closely, that risked breaking the "must
be designed from scratch" requirement. We threw it out and rebuilt it using o— and the replacement worked better.

We checked whether predictions were reproducible. Running the structure predindom seeds, only about a quarter of designs put the binder in the same placeboth times. Those that didn't were dropped. This is the check that most changed the final list.

---

The two objectives

Binding — this is the solid part. Our best designs score as well as the natural hormone does, and clearly better than control proteins that have no reason to bind. We also checked the
predicted poses are physically sensible: the binders sit in the water-filled buried in the greasy membrane where a soluble protein couldn't go.

Activation — this is the honest long shot. Activation is triggered when sometdeep in the receptor called E387. So we aimed designs at the deep pocket, andranked designs that reach E387 above those that don't. Nine of the 20 touch it.

We also found a clean reason why the other ten probably won't activate: the outer domain they bind doesn't change shape when the receptor switches on. So a protein gripping it has no
way to flip the switch — which explains why every published binder to that reot an activator.

---
                                                                                                                                                                                        mitations — the real ones
                                                                                                                                                                                        No design is predicted to activate. Touching E387 is a structural argument,edicts the cellular response.
- One model, not two. We only had Boltz-2 for the final scoring. Two runs of the same model agreeing means it's consistent, not that it's right. An independent second model would have been much stronger evidence.
- Missing software. No Rosetta licence, so we skipped a refinement step the standard pipeline normally includes. We measured the cost of that: it noticeably degrades the scores, which
  is why borrowed thresholds didn't transfer.
- Only 7 of 20 are fully reproducible. The other 13 are "partial" — the prediction shifts somewhat between runs.
- One thin control. Our "known good" benchmark rests on essentially a single

---

Pros and cons

In its favour
- Every design passed multiple independent filters, each calibrated on controls rather than assumption
- Two different target sites, so a wrong assumption about one doesn't sink ev
- Designs physically check out — right place, no membrane burial, no clashes, and the attachment tag used in the assay isn't blocked
- Fully documented: what we did, what we didn't, and what we got wrong along

Against it
- Activation remains a hope with a structural rationale, not a prediction
- Confidence rests on one model family
- Prediction scores are proxies — no actual binding strength was predicted for any design
- We had ~32 hours on one GPU; the published successes in this area used far