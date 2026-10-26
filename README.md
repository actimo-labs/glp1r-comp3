# GLP1R Track 3 — Anthropic × Adaptyv Challenge 03

Upload package, 2026-10-26. Layout mirrors the challenge-02 `tnf-upload` tree.

## The submission file

`results/submission.csv` (identical copy: `results/glp1r_designs.csv`)

20 designs, rank-ordered, top row highest. Columns are exactly those required by
FAQ 5: `name`, `sequence`, `molecule_class`.

Verified before packaging — all checks pass:

| Check | Result |
|---|---|
| ≤ 20 designs (Track 3 cap) | 20 |
| columns exactly name/sequence/molecule_class | yes |
| length 10–250 aa | 65–85 |
| single chain | yes |
| molecule_class | `single_chain` |
| 20 standard residues only | yes |
| de novo: no run > 5 vs GLP-1 / exendin-4 / glucagon / GIP / oxyntomodulin | longest run 5 |
| no HAEG motif | yes |
| no pair above 70% identity | max 65.3% |
| no Twin-Strep tag in the sequence | yes (organisers append it) |
| rank-ordered | row 1 = classAL85_2_s4 |
| target sequence | PDB 6X18 chain R; differs from UniProt P43220 only at L260F, which is 58 Å deep and >28 Å from every design |

`molecule_class` is optional and may also be set per design after upload. scFvs
and Fabs would be submitted as `{VH}:{VL}`; ours are all single-chain proteins,
so no sequence contains a colon.

## Why `molecule_class = single_chain`, verified

`molecule_class` describes **format**, not size. `nanobody` means a VHH domain —
a specific inherited immunoglobulin scaffold — not simply "a small
single-domain binder". The designs here were checked against that definition
rather than assumed:

| Evidence | These designs | A nanobody (VHH) |
|---|---|---|
| Secondary structure | **87% helix, 0% strand** (mean, Cα-based) | ~45–50% strand, <5% helix |
| Length | 65–85 aa | typically 110–130 aa |
| Cysteines | **0 in every design** | ≥2, for the invariant Ig disulfide |
| Framework motifs (`QVQLVESGGG`, `WYRQ`, `YYCA`, `WGQGT`) | **none present** | all characteristic |

The structural evidence is decisive: a nanobody is an immunoglobulin
β-sandwich, while these are α-helical bundles. Opposite folds.

Note also that ProteinMPNN was run with `--omit_AAs XC`, excluding cysteine by
construction — the canonical Ig disulfide was never buildable in this pipeline.

These are de novo RFdiffusion backbones with ProteinMPNN sequences and no
parental scaffold of any kind, which is consistent with the de novo requirement
and with the sequence check finding no run longer than 5 residues against any
natural protein tested.

## Composition

10 class A (TMD orthosteric vestibule) + 10 class C (ECD peptide-binding
groove), from 19 distinct backbones. **All 10 class A designs contact E387**,
the free-α-amine partner and activation trigger.

Scores are the mean of two Boltz-2 seeds. Designs are **tiered by pose
reproducibility**: ranks 1–10 place the binder within 3 Å across two
independent seeds, ranks 11–20 within 3–8 Å. Designs whose pose moved more than
8 Å between seeds were excluded outright (22 of 46 candidates).

Score agreement between seeds does **not** imply pose agreement — one candidate
agreed to within 0.03 on interface confidence while placing the binder 37 Å
apart. See `results/METHODS.md` §5a.

## Contents

| Path | What |
|---|---|
| `results/submission.csv` | **the file to upload** |
| `results/METHODS.md` | full methods, thresholds and provenance, not-computed list |
| `results/glp1r_metrics.csv` | per-design evidence |
| `results/rationale.md` | epitope, mechanism, evidence, failure mode per design |
| `results/final_ranked_full.csv` | the 20 with all metrics |
| `results/rejected_candidates.csv` | rejected designs with reasons |
| `results/selfconsistency.csv` | binder RMSD, prediction vs design — the gate that dropped class B |
| `results/boltz_seed{1,2}_scores.csv` | raw Boltz scores, both seeds |
| `results/af2_classC_scores.sc` | AF2-initial-guess scores, class C |
| `structures/` | 20 predicted complexes (binder chain A, receptor chain B) |
| `controls/` | positive and negative controls underpinning every threshold |
| `scripts/` | full pipeline |
| `target/` | prepared target, computed hotspots, epitope depth classes |
| `runlog/` | chronological record including incidents and corrections |

## Activation (objective 2) — where it is documented

| Question | Where |
|---|---|
| How activation entered the pipeline (epitope choice, target construction, selection, validation) | `results/METHODS.md` §4a |
| What the structural evidence supports, and what could not be computed | `results/METHODS.md` §9a, `runlog/08-activation-analysis.md` |
| Per-design measured trigger contacts | `results/rationale.md`, `results/glp1r_metrics.csv` |

In short: all 10 class A designs contact **E387**, the free-α-amine partner
whose engagement by GLP-1 His7 initiates activation, at a binder RMSD of
1.0–3.2 Å from the design (so this is the designed pose, not a repositioned
one). Their epitopes shift 3.8–8.6 Å between active and inactive receptor
states. The 10 class C designs target the ECD groove for affinity and are
**not** expected to activate.

## Honest summary

- **No KD is predicted.** Boltz confidence metrics are interface-plausibility
  proxies, not affinity estimates.
- **No design is predicted to activate.** Class A engages E387 and sits on a
  state-dependent epitope, which is the structural case for agonism — not
  evidence of it.
- **Class C is expected to bind and not activate.** The ECD is conformationally
  invariant between receptor states (1.21 Å internal RMSD vs 4.72 Å for the
  TMD), so an ECD binder has no route to shifting the equilibrium. This matches
  the published de novo GLP1R ECD binders being antagonists.
- No Rosetta (no licence), no MD, no RF3 cross-check. See `METHODS.md` §9 for
  the complete not-computed list.