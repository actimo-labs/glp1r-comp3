# Challenge 03 — GLP1R agonist: Claude Code handoff brief (v2)

Revised 2026-10-25. Supersedes v1. Source of truth: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/glp1r

**Deadline: Sun 25 Oct 23:59 AoE = Mon 26 Oct 22:59 Melbourne (AEDT). ~32 hours from the time this was written.**

Labels used throughout: **[OFFICIAL]** = competition rule; **[EVIDENCE]** = published experimental finding; **[HYPOTHESIS]** = structural proxy with no validation against experimental efficacy.

---

## 0. Working agreement for the agent

**Environment.** The session starts in `~/tnf-track3`. All GLP1R work goes in `~/tnf-track3/glp1r-comp3`. The tool checkouts and model resources (ProteinMPNN, RFdiffusion, Boltz, weights, envs) live alongside under `~/tnf-track3`.

**First action, before anything else:** `ls ~/tnf-track3`, then `git -C ~/tnf-track3/glp1r-comp3 rev-parse --show-toplevel` and `nvidia-smi`. Record in the run log: the actual directory layout, which tool checkouts exist and at what commit, which model weights are already on disk, the conda/venv in use, and the GPU state. Do not assume any path, checkout or checkpoint from this brief — verify and write down what is actually there. If a needed model or weight set is missing, say so immediately rather than downloading a large checkpoint unprompted.

**Git.** If `glp1r-comp3` is its own repository, work on `main` there. If it is a subdirectory of the `tnf-track3` repository, create and work on branch `glp1r-challenge03`. Either way: never commit onto the TNF-α working branch, and never modify or delete challenge-02 artifacts.

- Fresh session. The challenge-02 TNF-α work is a source of reusable *code* (manifest schema, QC, ranking, env), not a template for the *strategy* — the problem is structurally different.
- **Commit constantly and atomically.** Harjas has explicitly authorised commits. Commit after every meaningful unit: target prep done, a generation batch finished, a filter implemented, a metric computed, a shortlist cut. Message format: `glp1r: <what changed>` plus, where relevant, the counts (`glp1r: RFdiffusion batch 2, 480 backbones, seed 7`).
- Commit the code, configs, manifests, metrics tables, and run logs. Do **not** commit model weights, large PDB trajectory dumps, or raw structure sets over ~50 MB — put those under `.gitignore` and record their paths in the manifest instead.
- **No pushes to remote, no force-push, no history rewriting, no submission to the competition site without explicit approval from Harjas.** Commits are pre-authorised; everything that leaves the machine is not.
- Checkpoint every long batch so it resumes rather than restarts. With ~32 hours there is no budget for a lost run.

---

## 1. Verified challenge spec [OFFICIAL]

| Item | Value |
|---|---|
| Target | Human GLP1R, UniProt P43220-1, residues 24–432 |
| Assay construct | Sino Biological 13944-H18H-DD, DDM/CHS detergent, N-term FLAG + C-term His, **undisclosed thermostabilising mutations** |
| Structure reference | PDB 6X18 chain R |
| Recommended epitope | Peptide-binding groove of the ECD (hint, not a rule) |
| Length / format | 10–250 aa; single-chain protein, nanobody, scFv or Fab; de novo only |
| Cap | 20 designs (Track 3) |
| Binding gate | Design **immobilised** (C-terminal Twin-Strep), receptor in solution, pH 7.4. Counts as binding if the sensorgram fits a KD, or association shift > 300% over negative control |
| Activation | Only **confirmed binders** proceed. Cell-based CRE-SEAP reporter, Gαs → cAMP → luminescence, vs untreated and a reference agonist |
| Ranking | 1) affinity (lowest KD among all binders, groove-binding counts in its favour); 2) activation relative to reference agonist. Winner must do both |
| Selection for testing | Claude selects Track 2/3 designs from all submissions using a published prompt, with human approval. **Submitted metrics, structures and method descriptions are part of what it sees** |

Still to close before submitting: FAQ items 3–20 (use the page's Copy Markdown button — plain fetch truncates), selection-prompt tabs 02–09 (tab-gated, needs a browser), and the submission form's exact fields.

---

## 2. The scientific core: what actually makes an agonist

**Two-domain mechanism [EVIDENCE].** The ECD captures the ligand's C-terminal helix — that's the affinity, the "address". The ligand's N-terminus inserts into the 7TM core — that's the efficacy, the "message". Truncating just His7–Ala8 from GLP-1 costs ~100–300× affinity but ~1,000–10,000× cAMP potency. The free α-amine of His7 contacts E387(7.42b); the deep polar network includes M231(3.36b), Q234(3.37b), N240(3.43b), H363(6.52b), Q394(7.49b); ECL2 around W297 is essential.

**The decisive precedent [EVIDENCE].** Muratspahić, Feldman, Kim et al., *Nature* 656:1044–1053 (2026), "De novo design of miniproteins targeting GPCRs":
- They targeted the **GLP1R ECD** and got clean binders — 5.3 nM and 27 nM, IC50 39 and 61 nM — that were **antagonists**. ECD binding does not activate.
- They got **agonists** on class A receptors (MRGPRX1 EC50 390 nM and 1 µM, optimised to 42 nM; NK1R, 20 agonists spanning 1–231 nM) by targeting a **large TM2–TM7 orthosteric epitope on active-state structures**, explicitly reasoning that active-state stabilisation alone would be sufficient for agonism.
- **They used no MD, no ensemble sampling, no normal-mode analysis.** Filters were static: RFdiffusion → ProteinMPNN → Rosetta FastRelax → AF2 initial-guess, accepting on pAE_interaction, pLDDT_binder, Rosetta ddG, SAP. Agonism was confirmed only in cells.

**Therefore: design against the ACTIVE state (6X18 chain R), targeting the orthosteric vestibule, not the ECD alone.**

**On JJ's "design against the inactive state" suggestion.** The instinct behind it — that agonism is about conformational equilibrium, not interface score — is right and is the most useful idea in that conversation. But the specific inversion is the less-supported mechanistic model. GPCR pharmacology (β2AR ¹⁹F-NMR, pressure-resolved DEER, smFRET) favours **conformational selection**: agonists preferentially bind and stabilise a pre-existing active-like conformation rather than levering an inactive receptor open. Designing against the inactive or relaxed-apo surface optimises complementarity to the wrong geometry — it is how you make an antagonist. Keep the equilibrium framing; apply it as a *filter* on active-state-designed binders, not as the design target.

---

## 3. Agonism filters: what survives scrutiny

JJ's proposed metric — sample ~1000 conformations of each binder–receptor complex, count the fraction near the active state, use that as a scalar agonism score — is biophysically well-motivated. Ligand efficacy genuinely does track active-state population: DEER on β2AR resolves inverse agonist / antagonist / partial / full agonist by the fraction of TM6-outward population [EVIDENCE]. Four things break it in practice here:

1. **BioEmu-1 cannot do this at all.** It samples single soluble protein chains; the authors state membrane environments and ligands are neither represented in the model nor the training. BioEmu-2 (Oct 2026) *does* support multi-chain complexes and ligands and runs on an 80 GB A100 — but its own model card says sampling frequency **must not** be interpreted as a calibrated Boltzmann probability or free energy, which is precisely what the metric requires. No GPCR activation validation exists for either.
2. **MSA dependence.** A de novo binder is an MSA-less engineered construct — the regime where every MSA-dependent ensemble method degrades.
3. **The G-protein confounder — the serious one.** In 6X18 the active state is held by Gs and Nb35. Class B1 activation requires G-protein engagement; agonist alone does not fully stabilise the active state (directly shown for β2AR: agonist needs Nb80). Strip Gs for simulation and the intracellular face relaxes toward inactive **regardless of how good the binder is** — every design scores "antagonist", and the metric has near-zero discriminative power. Retaining mini-Gs or an α5 peptide fixes the collapse but introduces the opposite bias: the receptor is held active, so most bound designs look active. There is no clean escape, which is the fundamental reason this approach is fragile for class B1.
4. **Timescale.** Class B1 activation is millisecond. One A100 gives ~30–100 ns/day for a membrane-embedded GPCR complex after CHARMM-GUI setup and equilibration. That is three to six orders of magnitude short. The studies that did resolve GLP1R activation used supervised/biased CVs or Markov state models over weeks of GPU time. Unbiased MD is not a ranker here.

**Also: global backbone RMSD is the wrong order parameter**, even setting the above aside — it is dominated by ECD and loop motion, not the activation switch. Use class B1-specific metrics instead [EVIDENCE]:

| Order parameter | Active-state value |
|---|---|
| TM6 outward movement (Cα K346) | ~18 Å vs inactive |
| TM2(2x46)–TM6(6x37) Cα distance | GPCRdb class B inactive threshold ≤ ~13 Å; active well above |
| TM3–TM6 distance (R3.50–6.32) | split at 13.188 Å; agonist-bound median ~17 Å |
| PxxG kink angle (P6.47b–G6.50b) | ~90–112° active; inactive-like class B1 ~145–165° |
| Central polar network (K/R2.60b–N3.43b–H6.52b–Q7.49b) | intact |

Restrict all structural comparisons to the **TMD core**, excluding ECD and loops.

### The filter stack to actually build

**Stage 0 — design against active state.** 6X18 chain R, GLP-1 removed, active TMD geometry preserved, optionally mini-Gs/α5 retained to hold the intracellular face. RFdiffusion (motif-directed into the orthosteric vestibule) and/or BindCraft/BoltzGen → ProteinMPNN/SolubleMPNN → FastRelax → AF2 initial-guess. Baker's class-B thresholds as the starting point: pAE_interaction < 8–10, pLDDT_binder > 85–90, Rosetta ddG < −45, SAP < 60. **This stage has the strongest empirical support for producing agonists.**

**Stage 1 — per-state interface energetics [HYPOTHESIS].** PyRosetta FastRelax + InterfaceAnalyzer on each design modelled against active (6X18) and inactive (TMD core of 6LN2, or cleaner: the isolated TMD+NAM structures 5VEW/5VEX) receptors. Rank on ΔΔG_bind(active) − ΔΔG_bind(inactive). Minutes per design. Establishes relative complementarity per state; establishes nothing about dynamics or efficacy.

**Stage 2 — activation order parameters [HYPOTHESIS].** On the active-state complex, compute the table above with MDAnalysis/ProDy. Flag and demote any design that sterically clashes with the TM6-out active conformation — those are likely antagonists.

**Stage 3 — ANM/NMA tie-breaker [HYPOTHESIS].** ProDy ANM on apo-active, apo-inactive, and binder-bound receptor with the binder as added nodes. Does the binder keep the active state accessible along the low-frequency activation mode? Seconds per design. This is JJ's idea in its cheapest defensible form.

**Stage 4 — orthogonal structure prediction.** Boltz-2 and RF3, multiple seeds, agreement required on binding mode. Score on **inter-chain PAE and peptide/binder pLDDT, not ipTM** — on a 113-complex GPCR benchmark, median DockQ was AF2-ig 0.03, RF3 0.41, Boltz-2 0.56, and confidence metrics systematically overestimated correctness for misplaced peptides, worst in Boltz-2 where medium and high quality were indistinguishable by PAE alone. Note the opposite-direction caveat: AF2 is *active*-biased for class B1, so "AF2 predicts active" is weak evidence.

**Stage 5 — optional, only if time remains.** BioEmu-2 on the top handful, contact-conditioned toward active TM6-out geometry, read **qualitatively**. Report as uncalibrated. Do not present a population fraction as a probability of agonism.

**Nothing in stages 1–5 is a validated predictor of cellular cAMP activation.** The only validated readout is the assay. Say so in the methods write-up — the selection prompt rewards calibrated claims, and overclaiming an "agonism score" is the easiest way to look unserious.

---

## 4. Portfolio (20 designs)

Rebalanced from v1 toward pocket engagement, because the Baker GLP1R result shows ECD-only designs are antagonists.

| Class | n | Length | Rationale |
|---|---|---|---|
| **A. Orthosteric vestibule binders** (TM2–TM7 epitope, active-state designed) | 10 | 40–90 aa | The only architecture with a published route to agonism |
| **B. Two-domain mimetics** (N-terminal message module + ECD anchor) | 6 | 30–70 aa | Mirrors native mechanism; best affinity–efficacy balance if it folds |
| **C. ECD-groove minibinders** | 4 | 50–90 aa | Insurance for the binding gate and the affinity ranking; expect non-agonists |

Every design: activating element at the **free N-terminus**, Twin-Strep at the C-terminus. Pairs the native topology and the assay geometry. Any design whose mechanism needs a free C-terminus is dead on arrival.

No two submitted designs above 70% identity. Classes A and B must span at least three backbone topologies.

Antibody formats: not recommended. CDRs engage surface epitopes, not deep pockets; published GLP1R antibodies and Fabs are antagonists or neutral binders.

---

## 5. Open question that changes the portfolio — needs an answer before Stage 0 scales

The evidence says the reliable route to activation is reproducing GLP-1's N-terminal message (free α-amine into E387, His-Ala-Glu-Gly). The successful ProteinMPNN GLP-1RA campaign (Wei et al., *Adv. Sci.* 2025) fixed 13 conserved GLP-1 positions and redesigned 18, reaching EC50 ~0.011–0.036 nM.

**[OFFICIAL] says de novo designs only.** Grafting GLP-1's N-terminal residues onto a designed scaffold is arguably a violation on a strict reading, and arguably fine on a looser one — the scaffold is genuinely designed, and the motif is a pharmacophore rather than a known binder being mutated. This is a rules question, not a science question.

Harjas to ask in the competition Slack. Until it is answered, class B should be built so the message module can be swapped between a GLP-1-derived motif and a fully de novo one without regenerating the scaffolds.

---

## 6. 32-hour schedule

| Window | Work |
|---|---|
| 0–3 h | Layout/env/GPU audit in `~/tnf-track3`, branch setup in `glp1r-comp3`, target prep (6X18 chain R, two target definitions, hotspots derived from GLP-1 contacts in the structure — compute them, do not hardcode from literature), 6LN2/5VEW inactive anchors aligned on TMD core |
| 3–6 h | Pilot: ~10 designs end to end through every stage. Record wall-clock per design per tool. Do not scale before this passes |
| 6–20 h | Main generation + ProteinMPNN + AF2-ig triage |
| 20–28 h | Boltz-2/RF3 orthogonal validation + Rosetta per-state ddG + order parameters on the top ~100–200 |
| 28–31 h | Ranking, shortlist of 20, submission package, methods report, pre-submission audit |
| 31–32 h | Buffer |

If class A is slower than ~20 min/design end to end, cut sample counts, not validation.

---

## 7. Deliverables

- `submission/glp1r_designs.csv` — ranked 20, fields matching the form exactly (confirm from the form first).
- Per design: epitope targeted, intended mechanism, predicted KD evidence, main failure mode.
- Methods report: generation stack and versions, hotspot derivation, filter thresholds and their provenance, de novo compliance check, ranking rubric, and an explicit **"not computed / not run"** list.
- Predicted complex structures for the shortlist if the form accepts attachments.
- Rejected-candidate log with reasons.
- Pre-submission audit for Harjas. Nothing is uploaded without his approval.

---

## 8. Honest expectation

Binding is tractable — BindCraft averages ~46% experimental success, BoltzGen reports 66% on novel targets, and GLP1R ECD binders at single-digit nM already exist in the literature. Expect several confirmed binders among 20.

Agonism is not. De novo class B1 agonism has a handful of published successes, each with either a retained native message or a large screening funnel behind it. A 32-hour, single-GPU, 20-design campaign realistically targets **0–2 designs that both bind and activate**. Design for that honestly: no design in this submission is an agonist until the reporter assay says so, and the methods report should state exactly that.