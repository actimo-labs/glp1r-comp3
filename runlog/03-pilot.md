# Run log 03 — Pilot, measured per-design timings

All numbers measured on the A100-80GB this session unless marked *estimate*.

## Stage timings

| Stage | Input size | Measured |
|---|---|---|
| RFdiffusion class C, full ECD (101 res target) | 8 designs / 374 s | **46 s/design** |
| RFdiffusion class C, 15 Å crop (81 res) | 8 designs / 374 s | **46 s/design** |
| RFdiffusion class A, full TMD (283 res) | 8 designs / 788 s | **98 s/design** |
| RFdiffusion class A, 15 Å crop (150 res) | — | **~56 s/design** |
| ProteinMPNN (soluble, 4 seq/backbone) | 16 backbones → 64 seqs / 76 s | **1.19 s/design** |
| Sequence threading | 64 designs | 0.1 s total |
| AF2-ig, first design | 414-res complex | 121.9 s (one-off JAX compile) |
| AF2-ig, steady state | 414-res complex | **13.6 s/design** (11.4 s model) |

**Cropping helps class A (1.75×) and does nothing for class C.** The ECD was
already small — 101 → 81 residues — so there was nothing to win. comp-02
measured the same ratio on its target (1.84 → 1.09 min/backbone).

## AF2-ig correctness check

Positive control: native GLP-1 (chain P) + receptor (chain R) from 6X18, run
through the PyRosetta-free fork.

| Metric | Value |
|---|---|
| binder_aligned_rmsd | **0.496 Å** |
| target_aligned_rmsd | 0.791 Å |
| pae_interaction | 7.06 |
| plddt_binder | 84.0 |
| plddt_target | 93.4 |

Reproducing the cryo-EM peptide pose to 0.5 Å says the fork's feature path is
right. Chainbreak detection fired at residue 30 (the P/R boundary) and at the
two real gaps in chain R, which is the correct behaviour.

Re-running the same input four times gives byte-identical scores — the pipeline
is deterministic at `PRNGKey(0)`, so repeat sampling needs explicit seeds.

## Backbone yield — the real constraint

Interface QC (`scripts/12_interface_qc.py`) on the 16 uncropped pilot backbones.
Thresholds calibrated against the native complex: GLP-1, **backbone atoms only**
(the same representation RFdiffusion emits), contacts **28** receptor residues
at 5 Å. Full-atom it is 55, so measuring poly-GLY backbones roughly halves the
count and the threshold has to account for that.

| Class | ok | NOCONTACT | LIPID |
|---|---|---|---|
| C | 1/8 | 7 | 0 |
| A | 1/8 | 6 | 1 |

**~12% usable.** Low, but not anomalous — comp-02 screened ~500 backbones with
a 46% pass rate at its own calibrated cutoff, and its epitope was a flat
protein–protein interface rather than a membrane protein. One class A design did
land on the lipid-facing TM surface, which is the failure mode the filter exists
to catch, so the filter is earning its place.

## ProteinMPNN sanity check

First impression was alarming — sequences like
`...AAAALEERRRRELEARRAAAKAAGEAAVAAAAAVG...` look like low-complexity junk.
Checked against comp-02's known-good MPNN output on the same install:

```
comp-02 (worked):  SAEEHLRKARELVEEARAGVEAGSLSPEEITDLLMQAVSHLQQAAELLGA...   score 0.86
ours classC_1:     NEKAEKIAENLREALELAERAEELRRRAEELADDPELAARLEREAEELER...   score 1.03
ours classA_0:     MSMDAIDSFTNLAAAALEERRRRELEARRAAAKAAGEAAVAAAAAVGPAA...   score 1.10
```

Same E/A/R/L-rich character — that is simply what ProteinMPNN produces on de
novo helical bundles, not a bug. **MPNN score tracks backbone quality** and is a
usable cheap filter: comp-02's good designs sat at 0.86–0.93, and our worst
backbones are the ones producing poly-Ala runs at 1.10–1.23.

Per-backbone best score, pilot: classC_7 0.852 … classA_6 1.174. Note classC_7
has the best sequence score but **fails** interface QC, so the two filters are
complementary and both are needed.

## Bug found and fixed in my own code

`scripts/20_denovo_check.py` normalised pairwise identity by the aligned
overlap. A one-residue overhang that happened to match therefore scored 100%,
and every pair of designs looked like a duplicate — 560 spurious violations.
Now normalised by the shorter sequence with a minimum 80% overlap, and
self-tested in both directions (identical → 100%, 1-residue overhang → 0%,
contained-in-2× → 100%, 9/10 → 90%).

After the fix the only remaining identity failures are *sibling* sequences from
the same backbone (~72%), which is expected for 4 MPNN samples of one scaffold
and is handled by taking at most one design per backbone at selection time.

De novo checker self-test: native GLP-1 correctly FAILS (30-residue contiguous
match, plus the HAEG motif), a random sequence PASSES.

## Incident — challenge-02 artifacts overwritten and restored

Running comp-02's `scripts/02_crop_target.py` in place wrote
`target/crop_contigs.txt` and `target/crop_report.json` into the **comp-02**
`target/` directory. Its `--out` flag only redirects the PDB; those two paths
are derived from the script's own location.

Both files restored byte-identical from `backup_batch1/target/` and
`tnf-a-comp2/tnf-upload/target/`, with original mtimes preserved.
`target/tnf_AB_crop.pdb` was never touched. A sweep of everything modified under
`~/tnf-track3` on 2026-10-25 outside `glp1r-comp3/` found only files in the
`dl_binder_design` clone created this session.

Fix: the script is vendored into `glp1r-comp3/scripts/` so its ROOT resolves
inside this repo. Rule going forward: never execute a comp-02 script in place.

## Still outstanding

- Boltz-2 (`scripts/13_boltz.py`) written, template CIF built, **not yet run**.
- Class B not yet sampled. Input built: cropped receptor (243 res) + GLP-1 7–15
  backbone in native frame, with **all 27** message-module contacts retained.
- No FastRelax between MPNN and AF2 (no PyRosetta). Upstream dl_binder_design
  does relax. Goes in the "not computed / not run" list.
