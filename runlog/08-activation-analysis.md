# Run log 08 — Activation analysis

The goal is agonism, so this asks what structural evidence actually exists for
it. The answer is more negative than hoped, and the negative part is the more
useful finding.

## What could not be computed, and why

Brief §3 lists five class B1 activation order parameters. Two are not available
on this structure:

- **TM3(3.50b)–TM6(6.32b) distance.** 6.32b (~343) falls inside the receptor's
  disordered 339–343 gap and is absent from 6X18. Separately, the class A
  DRY-motif arithmetic does not transfer: GLP1R's TM3 at the derived 3.50b
  position is Glu, not Arg (sequence 235–260 is `YCVAANYYWLLVEGVYLYTLLAFSVF`,
  containing no Arg at all in 240–257).
- **TM2(2.46b)–TM6(6.37b) distance.** 6.37b (348) is present, but 2.46b cannot
  be placed without GPCRdb numbering, and guessing it by arithmetic is how the
  3.50b error above would have been made silently.

Verified and available: **P6.47b = P358** and **G6.50b = G361** (the TM6 `PLLG`
PxxG motif), **K346** as the TM6 reference, and the polar-network residues
N240 / H363 / Q394.

## Why order parameters on our own predictions prove nothing

Boltz pins the receptor to the **active-state** template (6X18 chain R).
Every predicted complex therefore shows active geometry regardless of what the
binder does. Measuring TM6 position on those is circular — the state is
imposed, not discovered. This is the same confounder the brief identifies for
MD: strip Gs and the receptor relaxes toward inactive no matter how good the
binder is; retain it and everything looks active.

So the question was reframed to a non-circular, purely steric one: **does each
design's epitope even change between receptor states?** Under conformational
selection, an agonist must preferentially engage a conformation that is
distinct in the active state. If the epitope is identical in both, the design
has no structural route to shifting the equilibrium.

## Reference validation

Inactive reference **6LN2** (full-length GLP1R with a negative allosteric
modulator), superposed on 6X18 chain R over 277 shared TMD-core Cα
(RMSD 4.72 Å; ECD and loops excluded, per the brief's note that global RMSD is
dominated by ECD motion).

**TM6 K346 Cα displacement between states: 17.4 Å** — matching the ~18 Å the
brief cites for TM6 outward movement. The reference and the superposition are
sound.

## Result 1 — no design clashes with either state

All 20 submitted designs show **zero steric clash** against both the active and
the inactive receptor. The reason is geometric: TM6's outward swing is on the
**intracellular** face, while every design binds the **extracellular**
vestibule or ECD. Extracellular binders simply do not overlap the activation
switch, so steric state-selection is not available to them.

## Result 2 — the ECD is conformationally invariant

Internal RMSD of each domain between active and inactive, each fitted on itself:

| Domain | Internal RMSD |
|---|---|
| **ECD (29–129)** | **1.21 Å** |
| TMD (136–423) | 4.72 Å |

**The ECD is rigid.** Between receptor states it reorients as a body but does
not change shape. A binder on the ECD groove therefore sees an identical
surface in the active and inactive receptor and has **no structural basis to
shift the conformational equilibrium**.

This is a direct mechanistic explanation for the result the brief quotes — that
the published de novo GLP1R ECD binders, at 5.3 nM and 27 nM, were
**antagonists**. It is not bad luck or a design flaw; the epitope carries no
state information.

**Consequence for our submission: the 10 class C designs should be expected to
bind and not to activate.** They are included for the binding gate and the
affinity ranking, which is what §4 of the brief allocated them for.

## Result 3 — class A epitopes do change with state

Median Cα displacement of the receptor residues each design actually contacts:

| Class | epitope shift (median across designs) |
|---|---|
| C (ECD groove) | 13.4 – 20.5 Å — but this is ECD **rigid-body reorientation**, not a shape change (see Result 2) |
| A (TMD extracellular) | **3.8 – 8.6 Å**, genuine local conformational change |

Class A epitopes sit on the TMD, which has a real internal conformational
difference between states (4.72 Å). The shifts are modest, but they are the only
activation-relevant structural coupling anywhere in this submission.

Combined with the finding that 7 of 10 class A designs contact **E387** — the
free-α-amine partner and activation trigger — class A carries whatever agonist
prospect this submission has.

## Honest conclusion

- **No design is predicted to activate.** Nothing here forecasts cAMP output.
- The structural case for agonism rests on class A: they engage the deep pocket
  including E387, and their epitope is state-dependent.
- The structural case *against* agonism for class C is strong and specific: a
  rigid, state-invariant epitope.
- The whole analysis is a **structural proxy**. It is not calibrated against any
  experimental efficacy measurement, and no agonism score is claimed.
