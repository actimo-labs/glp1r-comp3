# Per-design rationale

Submitted sequences carry no Twin-Strep tag; the C-terminal capture tag is assay infrastructure added by the organisers. Every design places its functional element at the N-terminus and tolerates a C-terminal fusion.

## 1. classAL85_2_s4

- **length**: 85 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 5 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233, T391, Y241**. Binder RMSD to the design is 2.07 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 6.71 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 1.71 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 8.04, ipTM 0.899. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 2. classAL65_9_s3

- **length**: 65 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 4 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233, T391**. Binder RMSD to the design is 3.16 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 21.24 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 2.36 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 8.93, ipTM 0.854. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 3. classAL65_24_s1

- **length**: 65 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 2 activation-trigger residue(s) in the Boltz-2 complex: **Q234, M233**. Binder RMSD to the design is 2.39 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 15.19 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 2.8 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 8.95, ipTM 0.862. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 4. classAL85_37_s3

- **length**: 85 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: no activation-trigger residue contact detected; submitted on binding evidence alone.
- **immobilisation geometry**: C-terminus is 7.3 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 0.62 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 10.3, ipTM 0.801. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 5. classCL65_39_s3

- **length**: 65 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 10.43 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 1.4 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 12.02, ipTM 0.878. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 6. classCL85_47_s4

- **length**: 85 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 7.27 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 1.71 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 14.84, ipTM 0.911. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 7. classAL75_39_s4

- **length**: 75 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 5 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233, T391, Y241**. Binder RMSD to the design is 2.68 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 17.38 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 3.75 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 9.41, ipTM 0.904. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 8. classAL75_19_s1

- **length**: 75 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 3 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233**. Binder RMSD to the design is 1.76 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 27.18 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 3.08 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 9.46, ipTM 0.888. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 9. classAL85_39_s4

- **length**: 85 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 3 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233**. Binder RMSD to the design is 2.41 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 20.66 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 5.06 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 10.37, ipTM 0.866. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 10. classAL85_39_s2

- **length**: 85 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 3 activation-trigger residue(s) in the Boltz-2 complex: **E387, M233, T391**. Binder RMSD to the design is 2.45 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 16.59 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 6.82 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 10.54, ipTM 0.825. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 11. classAL85_36_s3

- **length**: 85 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 5 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233, T391, Y241**. Binder RMSD to the design is 3.57 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 18.54 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 5.63 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 10.6, ipTM 0.892. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 12. classAL85_4_s3

- **length**: 85 aa
- **epitope targeted**: orthosteric vestibule / deep TMD pocket
- **intended mechanism**: Folded miniprotein targeting the GLP-1 message epitope on the transmembrane domain. RFdiffusion was conditioned on the vestibule rim, but the generation target retained the full deep pocket (receptor residues 187-241 and 361-397), and the designs reach further than their backbones suggest: modelled with side chains they engage the deep activation-trigger residues. **This is the architecture carrying the activation objective.** Intended mechanism: occupy the orthosteric pocket the native agonist uses and contact E387, the free-alpha-amine partner whose engagement by GLP-1 His7 initiates activation.
- **activation rationale (measured)**: contacts 3 activation-trigger residue(s) in the Boltz-2 complex: **E387, Q234, M233**. Binder RMSD to the design is 3.07 A, so this is the designed pose, not a repositioned one. The epitope it occupies shifts 3.8-8.6 A between active and inactive receptor states, so it is state-dependent.
- **immobilisation geometry**: C-terminus is 18.45 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 7.38 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 11.42, ipTM 0.779. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Agonism is unverified. Contacting E387 is not the same as presenting a free alpha-amine to it in the native geometry, and no structural method here predicts cAMP output.

## 13. classCL65_12_s1

- **length**: 65 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 6.44 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 7.79 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 11.78, ipTM 0.853. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 14. classCL65_23_s1

- **length**: 65 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 18.9 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 5.54 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 14.27, ipTM 0.818. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 15. classCL85_13_s4

- **length**: 85 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 23.53 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 4.0 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 14.44, ipTM 0.895. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 16. classCL85_37_s1

- **length**: 85 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 9.05 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 5.29 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 15.58, ipTM 0.938. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 17. classCL75_2_s3

- **length**: 75 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 30.61 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 6.38 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 16.23, ipTM 0.749. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 18. classCL65_54_s1

- **length**: 65 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 16.94 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 6.71 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 16.49, ipTM 0.833. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 19. classCL65_47_s1

- **length**: 65 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation geometry**: C-terminus is 9.44 A from the receptor, pointing away from the interface, so the C-terminal Twin-Strep tag is free for chip capture.
- **pose reproducibility**: two independent Boltz seeds place the binder 4.45 A apart after superposing on the receptor (**partial**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 16.99, ipTM 0.733. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

## 20. classCL85_40_s2

- **length**: 85 aa
- **epitope targeted**: ECD peptide-binding groove
- **intended mechanism**: Miniprotein against the extracellular-domain groove that binds GLP-1's C-terminal address helix - the epitope the competition recommends. Intended mechanism: high-affinity capture of the receptor, targeting the affinity ranking and the binding gate.
- **activation rationale (measured)**: none. The ECD epitope is state-invariant (see failure mode).
- **immobilisation risk**: the C-terminus sits 3.08 A from the receptor in the predicted complex. The design is immobilised on the chip via a C-terminal Twin-Strep tag, so a buried C-terminus risks the tag being inaccessible or immobilisation sterically blocking the interaction. **Flagged; ranked last for this reason.**
- **pose reproducibility**: two independent Boltz seeds place the binder 1.78 A apart after superposing on the receptor (**yes**). Designs whose pose moved more than 8 A between seeds were excluded from the submission entirely. Note that score agreement between seeds does NOT imply pose agreement - this is measured separately.
- **predicted-KD evidence**: Boltz-2 pae_interaction 15.0, ipTM 0.62. No KD is predicted; these are interface-confidence proxies calibrated against folded non-binders, not affinity estimates.
- **main failure mode**: Expected to bind but NOT to activate. The ECD is conformationally invariant between receptor states (1.21 A internal RMSD, vs 4.72 A for the TMD): it reorients rigidly but does not change shape, so an ECD binder sees the same surface in the active and inactive receptor and has no structural route to shifting the equilibrium. This matches the published de novo GLP1R ECD binders (5.3 nM, 27 nM) being antagonists.

