# Run log 02 — Epitope depth analysis (and why it challenges the portfolio)

`scripts/01_epitope_depth.py`, 2026-10-25 ~04:50 UTC. Output:
`runlog/01_epitope_depth.log`, `target/epitope_classes.json`.

## Why I ran this

Two GPCR-specific failure modes sit between the computed hotspots and a usable
RFdiffusion input:

1. RFdiffusion has no concept of a membrane, so it will happily design binders
   onto the lipid-facing TM surface — a big hydrophobic convex patch, exactly
   what it likes. Those cannot bind in a detergent micelle.
2. A folded 60–90 aa scaffold cannot reach into a deep pocket. Passing deep
   pocket residues as `ppi.hotspot_res` asks for something geometrically
   impossible.

## Method

Membrane normal = first principal axis of TMD Cα (62% of variance; the loops
dilute it). Extracellular direction set by which side the ECD centroid falls on
— derived, not assumed (+52.2 Å).

Depth from the most extracellular TMD Cα turned out to be a bad ruler: that Cα
is an extracellular loop tip, not the vestibule mouth. It classified **all 27**
class A hotspots as unreachable, which is an artefact.

**GLP-1 itself is the honest ruler** — whatever depth the native ligand reaches,
a ligand can reach. Re-measured against it:

| Reference | Depth |
|---|---|
| GLP-1 Arg36 (C-term, in ECD groove) | −4.4 Å |
| GLP-1 Phe28 | 2.4 Å |
| GLP-1 Leu20 | 14.3 Å |
| **GLP-1 Asp15 — top of the message module, taken as the vestibule mouth** | **20.6 Å** |
| **GLP-1 His7 — deepest insertion, the activation trigger** | **31.2 Å** |

A receptor residue is "rim" if the shallowest GLP-1 residue touching it is at or
above the mouth; "deep" otherwise.

## Result

**Only 6 of 27 class A hotspots are rim: 141, 201, 298, 299, 380, 384.**
The other 21 are deep, including the ones that matter most:

| Residue | Depth | Role |
|---|---|---|
| **E387** | **35.1 Å** | the free-α-amine partner — the activation trigger |
| Q234 | 29.3 Å | direct peptide contact |
| M233 | 29.0 Å | deep pocket |
| Y241 | 39.5 Å | deepest |

384 is borderline — classified rim only because GLP-1 Asp15 sits exactly at the
mouth and reaches it at 3.91 Å.

No class A hotspot is lipid-facing, so that failure mode did not bite here. For
later filtering: **76 of 283 TMD residues are lipid-facing overall** and any
design whose interface lands there should be rejected.

## What this means — flagging, not improvising

The GLP-1 message epitope is a deep narrow cleft, not a surface. A folded
miniprotein can engage a 6-residue rim patch and **cannot reach E387 at 35 Å**.
That is a problem for brief §4's portfolio, which puts **10 of 20 slots** on
class A ("the only architecture with a published route to agonism").

The published route it cites is *class A receptors* — MRGPRX1, NK1R — whose
orthosteric pockets are shallow and wide. Class B1 GLP1R's is neither. The same
paper targeted the GLP1R **ECD** and got clean binders that were **antagonists**.
So neither of the Baker results supports a folded miniprotein in the GLP1R
vestibule; that architecture is untested, and this measurement says it cannot
reach the trigger.

The architecture that *can* reach 35 Å is an extended N-terminal segment on an
anchored scaffold — the native topology, i.e. brief §4 **class B**, currently
allocated 6 slots.

Reading of the three classes on this evidence:

| Class | Can it bind? | Can it activate? |
|---|---|---|
| **C — ECD groove** (17 hotspots, depth ≤ 0 Å, fully exposed) | Yes — published 5.3 / 27 nM | No. Published ECD binders are antagonists |
| **B — two-domain mimetic** | Plausible via the ECD anchor | **Only architecture that reaches E387** |
| **A — vestibule-rim miniprotein** | 6-residue epitope is thin for high affinity | Does not reach the trigger |

Recommendation is to rebalance away from class A. Raised with Harjas rather
than actioned unilaterally — it is a departure from the brief.
