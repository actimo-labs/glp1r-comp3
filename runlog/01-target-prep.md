# Run log 01 — Stage 0 target prep

`scripts/00_prep_target.py`, run 2026-10-25 ~04:45 UTC. Full output in
`runlog/00_prep_target.log`, machine-readable result in `target/hotspots.json`.

## Source

`target/6X18.pdb`, fetched from RCSB this session. Cryo-EM, 2.10 Å, active
state. Chains A/B/G = Gs heterotrimer, N = Nb35, P = GLP-1(7-36)NH2,
R = receptor. Matches the brief exactly.

Chain R models **384 residues, 29–423**, with 11 internal gaps. The challenge
spec says the assay construct is UniProt P43220-1 **24–432**, so the structure
is missing 24–28 and 424–432 plus the gaps. Those are termini and loops, not
epitope, but the submission's method description should state the modelled range
rather than implying 24–432 was used.

## Numbering is clean — verified, not assumed

Before trusting any hotspot I checked chain R for a numbering offset against the
Wootten class-B1 positions quoted in brief §2, by asking what residue actually
sits at each number:

| Literature | Expected | At that number in 6X18 R |
|---|---|---|
| Q234 (3.37b) | Q | **Q** ✓ |
| N240 (3.43b) | N | **N** ✓ |
| W297 (ECL2) | W | **W** ✓ |
| H363 (6.52b) | H | **H** ✓ |
| E387 (7.42b) | E | **E** ✓ |
| Q394 (7.49b) | Q | **Q** ✓ |
| M231 (3.36b) | M | L — **Met is at 233**, see below |

So 6X18 chain R uses UniProt P43220 numbering and no offset correction is needed.

## Two findings

**1. The brief's "M231(3.36b)" is off by two.** Residue 231 in 6X18 chain R is
Leu; the Met is at **233**, which the contact calculation independently picks up
as a message-module contact. Every other position in brief §2 checks out. Using
233 going forward. Low impact — 233 was already in the computed hotspot set — but
worth correcting in the methods write-up rather than repeating the brief's number.

**2. The deep polar network residues are mostly *not* GLP-1 contacts.** N240,
H363, Q394 and M233 show no heavy-atom contact with GLP-1 within 4.5 Å. This is
not an error: brief §3's own table describes K/R2.60b–N3.43b–H6.52b–Q7.49b as the
*central polar network*, i.e. an intramolecular activation switch, not a ligand
contact set. Only **E387** (the His7 α-amine partner) and **Q234** are direct
peptide contacts among the §2 list. The practical consequence is that these
residues should be used as an **activation order parameter** (Stage 2 — is the
network intact?), **not** as RFdiffusion hotspots. Targeting them as hotspots
would be pulling the binder toward residues the native agonist does not touch.

W297 likewise has no contact at 4.5 Å, but 298/299/300 do — ECL2 engages through
the adjacent stretch.

## Domain split — computed

The ECD/TMD boundary was found by scanning every contiguous split of the Cα
contact graph (contacts < 10 Å, sequence separation > 4) and minimising
cross-boundary contacts normalised by the smaller side. No sequence numbering
assumed.

**Boundary = residue 136** (normalised cross-contacts 0.28).

| Domain | Range | Residues |
|---|---|---|
| ECD | 29–129 | 101 |
| TMD | 136–423 | 283 |

This is consistent with the conventional GLP1R ECD (~24–145) without having
assumed it.

## Hotspots — derived from GLP-1 contacts

GLP-1 partitioned by the two-domain mechanism: **message** = 7–15 (inserts into
the 7TM core), **address** = 16–36 (ECD groove). Receptor residues within 4.5 Å
heavy-atom of each half:

- message module contacts 27 receptor residues — **27/27 fall in the TMD**
- address helix contacts 29 — **17/29 fall in the ECD** (the other 12 are at the
  TMD extracellular face, which is expected: the helix spans both)

The partition is clean, which is a good sign the domain split is right.

**Class A target** (`target/glp1r_TMD.pdb`, orthosteric vestibule), 27 hotspots:
```
141 144 145 148 152 190 194 197 201 230 233 234 237 241
298 299 300 306 309 310 313 372 380 384 387 388 391
```

**Class C target** (`target/glp1r_ECD.pdb`, groove), 17 hotspots:
```
30 31 32 35 36 39 43 67 68 69 88 89 90 91 121 123 128
```

E387 is present in the class A set, as it must be — it is the free-α-amine
partner that brief §5 makes the core of the de novo message module.

## Outputs

| File | Atoms |
|---|---|
| `target/glp1r_R_full.pdb` | 3138 |
| `target/glp1r_TMD.pdb` | 2342 |
| `target/glp1r_ECD.pdb` | 796 |
| `target/hotspots.json` | contacts + both target definitions |

## Not yet done

- Gs/Nb35 retention variant of the class A target (brief Stage 0 says
  "optionally mini-Gs/α5 retained to hold the intracellular face"). Not needed
  for RFdiffusion into the extracellular vestibule; decide before Stage 2.
- 6LN2 / 5VEW inactive anchors aligned on the TMD core (brief §6, 0–3 h window).
  Needed for Stage 2 order parameters, not for generation. Pending.
