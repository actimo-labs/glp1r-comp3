# GLP1R Track 3 — Methods

Anthropic × Adaptyv Challenge 03. Single A100-80GB.

Numbers below are **measured on this machine** unless explicitly marked as an
estimate. Where a published threshold was not reproducible here, that is stated
rather than papered over.

---

## 1. Generation stack and versions

| Component | Version / provenance |
|---|---|
| RFdiffusion | Docker `rosettacommons/rfdiffusion:latest`, image `33d504109509` |
| RFdiffusion weights | `Complex_base_ckpt.pt`, shipped in the image (2023-03-30) |
| ProteinMPNN |  (2023-06-27), **soluble** weights `v_48_020` |
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

Only **6 of 27** message-epitope residues sit at or above the vestibule mouth,
so a folded miniprotein's **backbone** cannot reach E387. That conclusion drove
the original design plan — and it turned out to be too pessimistic by roughly
one side-chain length. Measured with side chains modelled, the submitted class A
designs do contact E387. See §5a and §4a; this paragraph is retained because it
is what the plan was based on at the time.

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

## 4a. How the activation objective was incorporated

The competition scores two things: affinity among all binders, and activation
relative to a reference agonist, with the winner required to do both. Activation
entered the pipeline at four points.

**1. Epitope choice (design time).** Activation of a class B1 GPCR is initiated
in the 7TM orthosteric pocket, not on the ECD. The deep pocket residues were
therefore computed from the native complex (§3) and used as RFdiffusion
hotspots for the TMD-targeting classes: **E387** (the free-α-amine partner
engaged by GLP-1 His7), Q234, M233, T391, Y241. Class C deliberately targets the
ECD groove instead and is submitted for affinity, not activation.

**2. Target construction (design time).** The class A generation target was
cropped to a 15 Å shell around its hotspots, which **retained the full deep
pocket** (receptor residues 187–241 and 361–397). The pocket was therefore
present and accessible during diffusion even though the hotspot list named only
the rim — which is why these designs reach deeper than the original depth
analysis predicted.

**3. Selection (ranking time).** Within class A, designs contacting **E387**
are ranked ahead of designs with better interface scores that do not. All 10
submitted class A designs contact E387; most also contact Q234, M233, T391 and
Y241. Per-design contacts are listed in `rationale.md` and `glp1r_metrics.csv`.

**4. Validation (post hoc).** A state-selectivity analysis (§9a) asks whether
each design's epitope differs between the active and inactive receptor, which is
the structural precondition for shifting the conformational equilibrium under a
conformational-selection model. Class A epitopes shift 3.8–8.6 Å between states;
the class C epitope does not change shape at all.

**What this does not amount to.** None of the above predicts cAMP output. There
is no agonism score in this submission, and none should be inferred from the
rankings. The honest claim is narrow: class A designs occupy the pocket the
native agonist uses, contact the residue whose engagement initiates activation,
and do so on an epitope that is state-dependent. Whether that produces signalling
is for the reporter assay.

### Molecule class

All 20 are submitted as `molecule_class = single_chain`. This was verified
against the definition of the alternatives rather than inferred from how they
were made:

| Evidence | These designs | A nanobody (VHH) |
|---|---|---|
| secondary structure (Cα-based, mean) | **87% helix, 0% strand** | ~45–50% strand, <5% helix |
| length | 65–85 aa | typically 110–130 aa |
| cysteines | **0 in every design** | ≥2 for the invariant Ig disulfide |
| framework motifs (`QVQLVESGGG`, `WYRQ`, `YYCA`, `WGQGT`) | **none** | all characteristic |

A nanobody is an immunoglobulin β-sandwich; these are α-helical bundles — the
opposite fold. ProteinMPNN was additionally run with `--omit_AAs XC`, so the
canonical Ig disulfide was not buildable by construction. No sequence contains
a colon, so none can be misparsed as an scFv/Fab `{VH}:{VL}` pair.

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

## 5a. Pose reproducibility across seeds — and why score stability is not enough

Boltz is stochastic. Running it twice gives two independent placements of the
binder, and comparing them is a direct test of whether a predicted pose is
reproducible at all.

**Score agreement does not imply pose agreement.** This was measured, not
assumed, and it was the most surprising result of the validation:

| Design | Δpae between seeds | binder pose RMSD between seeds |
|---|---|---|
| `classCL65_35_s4` | **0.03** | **37.1 Å** |
| `classAL75_33_s3` | 0.71 | 8.4 Å |
| `classAL85_2_s4` | 0.27 | 1.7 Å |

Two seeds can agree on interface confidence to within 0.03 while placing the
binder 37 Å apart. An earlier version of this submission used a
seed-**stability** gate on |Δpae|, which measured the wrong quantity and would
have ranked `classAL75_33_s3` first despite its pose not reproducing.

**Across the 46-design two-seed pool:**

| Pose RMSD between seeds | Count |
|---|---|
| < 3 Å (reproducible) | 10 |
| 3–8 Å (partial) | 14 |
| > 8 Å (divergent) | 22 |

Only 10 of 46 predicted poses reproduce. **The 22 divergent designs were
excluded from the submission outright.** The submitted 20 are ordered with all
10 reproducible designs first (ranks 1–10), then 10 partial (ranks 11–20), and
`pose_rmsd_seeds` / `pose_reproducible` are reported per design in
`glp1r_metrics.csv` and `rationale.md`.

This is the strongest statement available about positional correctness, and it
is still a statement about *reproducibility*, not about truth: two seeds of the
same model agreeing does not make the pose right. No orthogonal structure
predictor was available (§9).

---

## 5b. Physical plausibility of the predicted pose

Confidence and reproducibility say nothing about whether a pose is physically
possible in *this* assay. Three target-specific checks:

**Membrane penetration.** GLP1R sits in a DDM/CHS micelle and a soluble binder
cannot bury itself in the detergent phase. The membrane slab was estimated from
the **TMD Cα only** — including the ECD skews the principal axis and misplaces
the slab badly enough to report ECD binders as lipid-embedded, which is
geometrically impossible. On the correct slab:

| Class | atoms inside the slab | of those, in lipid (>8 Å from any receptor atom) |
|---|---|---|
| A | 54–94 | **0** |
| C | 0 | 0 |

Class A designs descend the **aqueous vestibule**, surrounded by receptor
protein the whole way — the same route GLP-1 takes. None occupies lipid. Class C
binds the ECD, entirely above the membrane.

**C-terminus accessibility — an assay-format constraint.** The design is
immobilised on the chip via a **C-terminal Twin-Strep tag**, with receptor in
solution. If the C-terminus is buried in the interface, immobilising it risks
the tag being inaccessible or sterically blocking the interaction being
measured. C-terminus distance to the receptor is therefore reported per design
and used as the **primary ranking tier**.

**19 of the 20 submitted designs have a free C-terminus** (6.4–30.6 Å from the
receptor). One, `classCL85_40_s2` at 3.1 Å, is retained but ranked **last** and
flagged in `glp1r_metrics.csv` and `rationale.md`. Earlier rankings contained
five such designs, including at ranks 4 and 5.

**Interpenetration.** Severe clashes (<2.2 Å) are 0–5 heavy-atom pairs across
all 20 — consistent with normal close packing, no interpenetration.

---

## 5c. Self-consistency — the gate that decided the submission

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

0. **C-terminus accessible** — ≥5 Å from the receptor, so the Twin-Strep tag
   used for chip immobilisation is free (§5b). Primary ranking tier.
0b. **pose reproducibility** — binder pose RMSD < 8 Å between two independent
   Boltz seeds. 22 of 46 candidates failed this and were excluded (§5a).
1. **self-consistency** — binder RMSD < 5 Å between prediction and design
2. **Boltz interface** — `pae_interaction` < 20.0 and `ipTM` > 0.56
3. **epitope** — interface not predominantly lipid-facing
4. **de novo** — no run > 5 residues matching a natural incretin peptide

Scores are the **mean of two Boltz seeds**, not one. Single-seed ranking proved
noisy: mean |Δpae| between seeds was 1.91, maximum 8.85.

Ranking is tiered by **pose reproducibility first** (§5a), then by score within
tier. Ranks 1–10 reproduce their pose across seeds (< 3 Å); ranks 11–20 are
partial (3–8 Å). Designs whose pose diverged by more than 8 Å were excluded.

Then ranked on mean `pae_interaction` ascending, `ipTM` descending as tie-break.
`plddt_binder` is **not** used: it was shown three separate times to track
foldedness rather than binding.

Within class A, designs contacting **E387** are preferred before ranking on
score, because activation is the objective and E387 is the free-α-amine partner.
All 10 submitted class A designs contact it.

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
  The 10 class A designs contact E387 and sit on a state-dependent epitope,
  which is a structural *case* for agonism, not evidence of it. Whether any
  design activates is unknown until the reporter assay says so.
- **No Rosetta FastRelax, ddG, or SAP** — no PyRosetta licence.
- **No per-state interface energetics** (planned Stage 1) — not run.
- **No MD, no ensemble sampling, no BioEmu** — not run. The published rationale
  for skipping it holds: class B1 activation is a millisecond process, an A100
  gives ~30–100 ns/day for a membrane GPCR complex, and stripping Gs collapses
  the active state regardless of binder quality.
- **No activation order parameters computed on our predicted complexes.** Not
  an oversight: Boltz pins the receptor to the active-state template, so every
  prediction shows active geometry regardless of the binder and the measurement
  would be circular (§9a). Two of the five planned parameters are also
  uncomputable on this structure — 6.32b (~343) lies in the disordered 339–343
  gap, and GLP1R's TM3 has no Arg at the derived 3.50b, so class A DRY-motif
  arithmetic does not transfer.
- **No RF3 cross-check.** Not installed. Stage 4 is Boltz-2 only, so
  "orthogonal agreement between two structure predictors" was not achieved.
- **No orthogonal structure predictor.** Pose reproducibility is measured
  between two seeds of the *same* model (Boltz-2). Agreement between seeds of
  one model is weaker evidence than agreement between two independent models,
  and RF3 was not installed (see below).
- **Target sequence.** Designs were made against **PDB 6X18 chain R**, not the
  UniProt sequence directly. 6X18 chain R differs from UniProt P43220 at
  **one position: L260F**, presumably a construct mutation for structure
  determination. Residue 260 sits at 58 Å depth on the intracellular side, is
  not a hotspot, does not contact native GLP-1, and **no submitted design comes
  within 28 Å of it**, so the difference cannot affect any designed interface.
  The assay construct's own thermostabilising mutations are undisclosed.
- **One valid positive control.** The exendin-4 positive failed (27.33) because
  it was threaded onto GLP-1's backbone rather than its own — a non-native
  pairing. The calibration ceiling rests on a single positive.
- **Designs span 29–423 of the receptor**, not the full 24–432 assay construct.

---

## 9a. Activation — what the structural evidence actually supports

Full analysis in `runlog/08-activation-analysis.md`. Summary:

**The ECD is conformationally invariant between receptor states.** Internal
RMSD active vs inactive (6LN2, superposed over 277 TMD-core Cα; TM6 K346
displacement 17.4 Å, matching the cited ~18 Å):

| Domain | internal RMSD |
|---|---|
| ECD (29–129) | **1.21 Å** |
| TMD (136–423) | 4.72 Å |

The ECD reorients as a rigid body but does not change shape. A binder on the ECD
groove therefore sees an identical surface in both states and has **no
structural route to shifting the conformational equilibrium**. This is a
mechanistic explanation for the published de novo GLP1R ECD binders (5.3 nM,
27 nM) being **antagonists** — the epitope carries no state information.

**Consequence: the 10 class C designs are expected to bind and not to activate.**
They are submitted for the binding gate and the affinity ranking.

**Class A carries the agonist case.** Their epitopes lie on the TMD, which does
differ between states (4.72 Å internal), and the specific residues they contact
shift 3.8–8.6 Å on activation. All 10 contact **E387**, most also Q234, M233,
T391 and Y241.

**What could not be computed.** Two of the five order parameters in the plan are
unavailable: 6.32b (~343) lies inside the receptor's disordered 339–343 gap, and
GLP1R's TM3 has no Arg at the derived 3.50b position, so class A DRY-motif
arithmetic does not transfer. Verified and usable: P358/G361 (PxxG), K346,
and the N240/H363/Q394 polar network.

Order parameters were **not** computed on our own predicted complexes, because
Boltz pins the receptor to the active-state template — every prediction looks
active regardless of the binder, so the measurement would be circular.

No design is predicted to activate. This is a structural proxy, uncalibrated
against any efficacy measurement.

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
