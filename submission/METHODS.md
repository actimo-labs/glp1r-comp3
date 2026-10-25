# GLP1R Track 3 — Methods

Anthropic × Adaptyv Challenge 03. Single A100-80GB, ~32 h wall clock.

Numbers below are **measured on this machine** unless explicitly marked as an
estimate. Where a published threshold was not reproducible here, that is stated
rather than papered over.

---

## 1. Generation stack and versions

| Component | Version / provenance |
|---|---|
| RFdiffusion | Docker `rosettacommons/rfdiffusion:latest`, image `33d504109509`, internal commit `fa340147b9006156b251d1ad0391e3ea8e5f73eb` (2025-06-22) |
| RFdiffusion weights | `Complex_base_ckpt.pt`, shipped in the image (2023-03-30) |
| ProteinMPNN | `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57` (2023-06-27), **soluble** weights `v_48_020` |
| AF2-initial-guess | `nrbennet/dl_binder_design` `cafa3853`, forked to remove PyRosetta (see §6) |
| AF2 weights | `params_model_1_ptm.npz` (2021-07-19) |
| Boltz-2 | pip `boltz` 2.2.1, checkpoints `boltz2_conf`, `boltz2_aff` |
| jax / jaxlib | 0.4.20 + CUDA 12.3-era NVIDIA libs (pinned; see §6) |

A host checkout of RFdiffusion at commit `86507b65` exists but is **not** what
executes — the Docker image is. The image commit is the one that matters.

---

## 2. Target preparation

Target: **PDB 6X18 chain R**, cryo-EM 2.10 Å, active state (Gs + Nb35 bound),
the structure reference named by the competition. Chain R models **384
residues, 29–423**, with 11 internal gaps. The assay construct is UniProt
P43220-1 residues 24–432, so 24–28 and 424–432 are not modelled; designs were
made against what the structure contains.

**Numbering was verified, not assumed.** Each class-B1 position quoted in the
literature was checked against the residue actually present at that number:
Q234, N240, W297, H363, E387, Q394 all match. One does not — the commonly cited
"M231 (3.36b)" is **Leu** in 6X18; the Met is at **233**. We use 233.

**Domain split was computed**, not taken from the literature. Every contiguous
split of the Cα contact graph (contacts < 10 Å, sequence separation > 4) was
scanned and the one minimising cross-boundary contacts, normalised by the
smaller side, was selected. **Boundary = residue 136.** ECD 29–129 (101 res),
TMD 136–423 (283 res). This agrees with the conventional ECD (~24–145) without
having assumed it.

---

## 3. Hotspot derivation

Hotspots were **computed from the native complex**, not copied from the paper.
Receptor residues within 4.5 Å heavy-atom of GLP-1 (chain P), partitioned by the
two-domain mechanism — message module = GLP-1 7–15, address helix = 16–36.

- message module contacts **27** receptor residues, all 27 in the TMD
- address helix contacts **29**, of which 17 in the ECD

A finding that changed the design: the deep polar network residues (N240, H363,
Q394, M233) make **no** contact with GLP-1 at 4.5 Å. They are an *intramolecular*
activation switch, not a ligand contact set. Using them as RFdiffusion hotspots
would pull binders toward residues the native agonist never touches. They are
used as activation order parameters instead.

### Epitope depth — why the portfolio is weighted as it is

Depth was measured along the TM bundle principal axis (62% of Cα variance), with
the extracellular direction set by which side the ECD centroid falls on. The
honest ruler is GLP-1 itself:

| Reference | Depth |
|---|---|
| GLP-1 Arg36 (ECD groove) | −4.4 Å |
| GLP-1 Asp15 — vestibule mouth | 20.6 Å |
| GLP-1 His7 — deepest insertion | 31.2 Å |
| **E387 — the free-α-amine partner, the activation trigger** | **35.1 Å** |

Only **6 of 27** message-epitope residues sit at or above the vestibule mouth.
A folded miniprotein can engage that thin rim but cannot reach E387. This is why
class A is treated as a binding-only architecture and why the agonist effort
went into longer binders that reach the pocket.

76 of 283 TMD residues are lipid-facing; designs whose interface lands there are
rejected, since they cannot bind in a detergent micelle.

---

## 4. Design classes

| Class | Epitope | Length | Generated | **Submitted** |
|---|---|---|---|---|
| A | vestibule rim (TMD) | 65/75/85 | 120 | **10** |
| B | deep orthosteric pocket | 90/110 | 250 | **0 — see §5a** |
| C | ECD peptide-binding groove | 65/75/85 | 198 | **10** |

The planned allocation was 4A / 10B / 6C. It did not survive validation: class B
failed self-consistency completely and was dropped, and class A turned out to be
the architecture that engages the activation trigger. Final split is 10A / 10C.

Binder lengths are **discrete, not sampled from a range**. JAX retraces AF2 for
every new complex length — measured 122 s on a new shape vs 10.5 s on a repeat —
so a continuous 60–90 range would have cost ~31 compiles per class. Three
discrete lengths cost 3. Length diversity is quantised, not removed.

Targets were cropped to a 15 Å shell around the hotspots for generation only
(class A 283→150 residues, 1.75× faster; class C 101→81, no gain). **All scoring
uses the complete 384-residue receptor** — designs are transformed back into the
original frame, which is exact because RFdiffusion's recentring is a pure
translation (measured: 214.278 Å shift, 0.133 Å RMSD after removing it).

---

## 5. Filters and thresholds, with provenance

**No published threshold is used unmodified.** Each was checked against controls
run through this exact pipeline, and two did not survive.

### Controls

- **Positive**: native GLP-1(7-36) — a validated sub-nanomolar binder and full
  agonist whose bound structure *is* the target.
- **Negative, unfolded**: composition- and length-matched shuffles.
- **Negative, folded**: six challenge-02 TNF binders — real de novo designed
  miniproteins, scored against GLP1R. FAQ 8 explicitly permits using known
  binders "to evaluate/calibrate design filtering metrics".

The folded negative is the one that matters, and it overturned a conclusion
drawn from the shuffle. A shuffled sequence does not fold, so it tests the wrong
failure mode.

### Interface QC

≥4 target residues within 5 Å of the binder (backbone-only representation), and
<35% of the interface on lipid-facing residues.

The threshold is calibrated against challenge-02's RFdiffusion backbones
measured identically: median 6.0, mean 6.3, max 13 — statistically the same as
this campaign's median 5.0. An earlier threshold of 8 would have rejected ~75%
of a backbone set that produced a real submission.

### AF2-initial-guess — class C only

`pae_interaction < 22.0`.

| Control | pae_interaction |
|---|---|
| native GLP-1, full atom | 7.06 |
| native GLP-1, backbone-only (our representation) | **17.39** |
| shuffled negatives (n=8) | min 22.00 |
| folded non-binders (n=3) | 26.10 – 27.45 |

**Baker's published `pae_interaction < 8–10` is not usable here.** In our
representation the native agonist scores 17.39, so that cutoff rejects a
validated full agonist. The gap is the missing Rosetta FastRelax step (§6).

**AF2-ig is not applied to classes A and B.** Their scores saturate at 27.3–28.4
— indistinguishable from, in fact marginally worse than, proteins known not to
bind. AF2-ig has median DockQ 0.03 on GPCR complexes; this is its documented
blind spot on membrane-protein interfaces, not evidence about the designs.

### Boltz-2 — the arbiter, all classes

Receptor pinned to the active-state template (6X18 chain R). Without a template
Boltz does not dock the chains reliably, and we are specifically designing
against the active state.

| Metric | POS native | class C best | folded non-binder | Used? |
|---|---|---|---|---|
| `pae_interaction` | 4.56 | 12.55 | **20.04 – 28.66** | **primary, < 20.0** |
| `ipTM` | 0.949 | 0.942 | **0.113 – 0.561** | secondary, > 0.56 |
| `plddt_binder` | 91.6 | 84.1 | 75.9 – 86.8 | fold check only, > 70 |

`plddt_binder` **overlaps the non-binder band entirely** and is not used for
ranking — it measures foldedness, not binding.

On `ipTM`: the brief and challenge-02's notes both advise against it, but both
concerned a *homotrimer* target where global ipTM is dominated by
target-internal chain pairs at ~0.9. GLP1R is a single chain, so ipTM *is* the
binder–target pair. Used as a secondary with that reasoning stated.

---

## 5a. Self-consistency — the gate that decided the submission

Boltz pins the receptor to a template but **re-predicts the binder pose**. A
good interface score therefore means either "binds as designed" or "Boltz found
some other pose it likes", and the brief's own benchmark warns that Boltz
confidence overestimates correctness for misplaced binders.

Separating those: superpose the predicted receptor onto the design's receptor,
then measure how far the binder moved.

| Class | binder RMSD vs design (median) | self-consistent (<5 Å) |
|---|---|---|
| C | **2.5 Å** | **52 / 60** |
| A | **3.1 Å** | **25 / 40** |
| B | **64.8 Å** | **0 / 100** |

**Class B was eliminated by this gate.** Its Boltz interface scores looked
survivable (7 of 100 under the PAE cut), but the predicted poses are ~65 Å from
the designed ones. Inspecting the best five: two are placed at 48–67 Å depth,
i.e. against the **intracellular** face, and the rest contact 59–76 receptor
residues spanning residues 29–392 — draped across the whole protein. Neither is
a credible binding mode. Submitting them would have meant 10 near-certain
non-binders selected on a confidence score that self-consistency contradicts.

A contributing weakness: the class B designed interfaces were thin to begin with
(6–23 contacts for a 90–110 aa binder). The interface-QC floor of ≥4 contacts was
calibrated on 65–85 aa binders and was too permissive for longer ones.

### Class A reaches the activation trigger — correcting §3

§3 concluded that a folded scaffold could not reach E387 at 35.1 Å. That was
measured on **backbone** geometry, and the RFdiffusion designs are backbone-only
poly-glycine. Modelled with side chains, which extend ~5 Å further, **7 of the
10 submitted class A designs contact E387**, and most also contact Q234, M233,
T391 and Y241 — at binder RMSD 1.0–3.2 Å from the design, so this is the
designed pose with side chains, not a different complex.

The depth analysis was too pessimistic by roughly one side-chain length. The
class A target crop retained residues 187–241 and 361–397, so the deep pocket
was present during generation even though the hotspots named only the rim.

This is the campaign's main scientific surprise: the architecture built *for*
the deep pocket (B) failed, and the one built for the rim (A) reached the
trigger. Stated plainly because it was not predicted.

---

## 6. Deviations forced by the environment

**No PyRosetta licence.** Consequences, all real:

1. `af2_initial_guess/predict.py` imports PyRosetta at module scope even for PDB
   input. Forked to remove it. The dependency was I/O only — the feature path
   round-trips pose → PDB text → coordinates, so reading the PDB directly gives
   identical features. Verified: the fork reproduces the native GLP-1 pose at
   **0.496 Å** binder RMSD.
2. **No FastRelax between ProteinMPNN and AF2.** Upstream relaxes and thereby
   builds side chains. Measured cost on the positive control: `pae_interaction`
   **7.06 → 17.39**, `plddt_binder` 84.0 → 51.5. This is why no published AF2
   threshold transfers.
3. Stage 1 of the planned filter stack (per-state Rosetta ddG via
   InterfaceAnalyzer) was **not run at all**.

A latent upstream bug was also fixed: the selenomethionine branch of
`af2_get_atom_positions` referenced an out-of-scope variable and would raise
`NameError` on any PDB containing MSE.

---

## 7. De novo compliance

FAQ 8: designs must be zero-shot, must not have "relied on an initial binder as
a starting point", and must have "adequate sequence- **and structural**-
diversity" from known proteins.

**A planned architecture was withdrawn on this rule.** Class B was first built
by holding GLP-1 residues 7–15 **backbone** fixed as an RFdiffusion motif, with
sequence designed de novo. That has full sequence diversity but deliberately
reduced *structural* diversity at the functional site, and GLP-1 is the native
binder of this exact target. It was discarded.

The submitted class B uses **no ligand coordinates of any kind**. RFdiffusion is
pointed at deep receptor residues (E387, Q234, M233, T391, Y241) and must find
its own way into the pocket. It performs better than the withdrawn version:

| | withdrawn motif version | submitted version |
|---|---|---|
| deepest reach | 30.4 Å | **43.7 / 44.4 Å** |
| activation-trigger residues contacted | 0 | **3 of 5** |
| E387 contact | none | **4.12 / 4.40 Å** |

Every submitted sequence is checked for contiguous runs >5 residues against
GLP-1(7-36), GLP-1(1-37), exendin-4, glucagon, GIP(1-42) and oxyntomodulin, and
for the HAEG motif. The check is duplicated inside the submission builder so it
cannot be bypassed.

Known binders are used only where FAQ 8 permits: as calibration controls.

---

## 8. Ranking rubric

Gates first, all of which a design must clear:

1. **self-consistency** — binder RMSD < 5 Å between prediction and design
2. **Boltz interface** — `pae_interaction` < 20.0 and `ipTM` > 0.56
3. **epitope** — interface not predominantly lipid-facing
4. **de novo** — no run > 5 residues matching a natural incretin peptide

Then ranked on `pae_interaction` ascending, `ipTM` descending as tie-break.
`plddt_binder` is **not** used: it was shown three separate times to track
foldedness rather than binding.

**Per-class quotas, 10 A and 10 C, rather than one global ranking.** Ranking
globally put 19 of 20 slots on class A. Two reasons to split instead: the
competition states that binding at the peptide groove (the ECD — class C)
counts in a design's favour, and a single-epitope submission fails entirely if
class A's predicted insertion turns out to be an artefact.

At most 2 sequences per backbone, no pair above 70% identity. The submitted 20
come from **20 distinct backbones**.

---

## 9. Not computed / not run

Stated plainly, because the selection prompt rewards calibrated claims.

- **No KD is predicted.** Boltz confidence metrics are interface-plausibility
  proxies, not affinity estimates. No design has a predicted binding constant.
- **No agonism prediction.** Nothing here forecasts cellular cAMP activation.
  Class B is *designed* to engage the activation trigger; whether it activates
  is unknown until the reporter assay says so.
- **No Rosetta FastRelax, ddG, or SAP** — no PyRosetta licence.
- **No per-state interface energetics** (planned Stage 1) — not run.
- **No MD, no ensemble sampling, no BioEmu** — not run. The published rationale
  for skipping it holds: class B1 activation is a millisecond process, an A100
  gives ~30–100 ns/day for a membrane GPCR complex, and stripping Gs collapses
  the active state regardless of binder quality.
- **No activation order parameters computed on final designs** (TM6 outward,
  TM2–TM6 distance, PxxG kink). MDAnalysis and ProDy are installed; this was
  deprioritised against validation.
- **No RF3 cross-check.** Not installed. Stage 4 is Boltz-2 only, so
  "orthogonal agreement between two structure predictors" was not achieved.
- **One valid positive control.** The exendin-4 positive failed (27.33) because
  it was threaded onto GLP-1's backbone rather than its own — a non-native
  pairing. The calibration ceiling rests on a single positive.
- **Designs span 29–423 of the receptor**, not the full 24–432 assay construct.

---

## 10. Honest expectation

Binding is plausible for class C: the ECD groove is the competition's
recommended epitope, single-digit nM binders to it exist in the literature, and
our best class C designs score at the native agonist's level against a folded
non-binder null.

Agonism is a long shot. De novo class B1 agonism has few published successes and
each had either a retained native message module or a large screening funnel.
This is a ~32 h, single-GPU, 20-design campaign without either. **No design here
is an agonist until the reporter assay says so.**
