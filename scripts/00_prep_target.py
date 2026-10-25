#!/usr/bin/env python3
"""
Stage 0 target prep for GLP1R challenge 03.

Inputs : target/6X18.pdb  (cryo-EM 2.10 A, active state: Gs + Nb35 + GLP-1 + receptor)
Outputs: target/glp1r_R_full.pdb   receptor chain R, cleaned
         target/glp1r_TMD.pdb      TMD core  (class A orthosteric vestibule target)
         target/glp1r_ECD.pdb      ECD       (class C groove target)
         target/hotspots.json      computed contacts + the two target definitions

Everything here is DERIVED FROM THE STRUCTURE. Per the brief, hotspots are
computed from the native GLP-1 (chain P) contacts, not copied from the paper.
The literature numbers are printed alongside only as a cross-check, and any
disagreement is reported rather than silently reconciled.

The two-domain mechanism gives us the two target definitions for free:
  GLP-1 His7..~Ala15   the "message" -> inserts into the 7TM core  -> TMD target
  GLP-1 ~Glu16..Arg36  the "address" -> ECD peptide-binding groove -> ECD target
So we partition GLP-1 by residue number, and each half's contact set defines
one epitope.
"""

import json
import os
import sys
from collections import defaultdict

import numpy as np

SRC = "target/6X18.pdb"
CONTACT_CUTOFF = 4.5   # A, heavy-atom
MESSAGE_MAX = 15       # GLP-1 7..15 = message module; 16..36 = address helix

# Cross-check only. Never used to define anything.
LIT_DEEP_POCKET = {231: "M3.36b", 234: "Q3.37b", 240: "N3.43b",
                   363: "H6.52b", 387: "E7.42b", 394: "Q7.49b"}
LIT_ECL2 = 297


def read_chain(pdb, chain):
    """Return {resnum: {'name': str, 'atoms': [(atomname, xyz)]}} for one chain."""
    res = {}
    for l in open(pdb):
        if l[:4] != "ATOM" or l[21] != chain:
            continue
        if l[16] not in (" ", "A"):
            continue
        elem = (l[76:78].strip() or l[12:16].strip()[:1])
        if elem == "H":
            continue
        n = int(l[22:26])
        r = res.setdefault(n, {"name": l[17:20].strip(), "atoms": []})
        r["atoms"].append((l[12:16].strip(),
                           np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return res


def contacts(a_res, b_res, cutoff):
    """{a_resnum: {b_resnum: min_dist}} for heavy-atom pairs under cutoff."""
    b_items = [(bn, np.array([xyz for _, xyz in b["atoms"]]))
               for bn, b in sorted(b_res.items())]
    out = defaultdict(dict)
    for an, a in sorted(a_res.items()):
        acoords = np.array([xyz for _, xyz in a["atoms"]])
        for bn, bcoords in b_items:
            d = np.linalg.norm(acoords[:, None, :] - bcoords[None, :, :], axis=-1).min()
            if d <= cutoff:
                out[an][bn] = round(float(d), 2)
    return out


def split_domains(res):
    """
    Find the ECD/TMD boundary from the Ca contact graph.

    The receptor is two structurally separable domains joined by a linker, so we
    scan every contiguous split point and keep the one that minimises contacts
    crossing the boundary (normalised by the smaller side, so it cannot just
    shave off a terminus). No sequence numbering is assumed.
    """
    nums = sorted(res)
    ca = {n: next((xyz for an, xyz in res[n]["atoms"] if an == "CA"), None) for n in nums}
    nums = [n for n in nums if ca[n] is not None]
    coords = np.array([ca[n] for n in nums])
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    adj = (d < 10.0)
    np.fill_diagonal(adj, False)
    # ignore near-neighbours in sequence: they contact trivially
    for off in range(1, 5):
        idx = np.arange(len(nums) - off)
        adj[idx, idx + off] = False
        adj[idx + off, idx] = False

    best, best_score = None, np.inf
    lo, hi = int(0.15 * len(nums)), int(0.85 * len(nums))
    for i in range(lo, hi):
        cross = adj[:i, i:].sum()
        score = cross / min(i, len(nums) - i)
        if score < best_score:
            best_score, best = score, i
    return nums[best], best_score


def write_pdb(src, chain, resnums, out, newchain="A"):
    keep = set(resnums)
    n = 0
    with open(out, "w") as f:
        for l in open(src):
            if l[:4] != "ATOM" or l[21] != chain:
                continue
            if l[16] not in (" ", "A"):
                continue
            if int(l[22:26]) not in keep:
                continue
            f.write(l[:21] + newchain + l[22:])
            n += 1
        f.write("TER\nEND\n")
    return n


def main():
    if not os.path.isfile(SRC):
        sys.exit(f"missing {SRC}")

    R = read_chain(SRC, "R")   # receptor
    P = read_chain(SRC, "P")   # GLP-1(7-36)
    print(f"chain R: {len(R)} residues  {min(R)}-{max(R)}")
    print(f"chain P: {len(P)} residues  {min(P)}-{max(P)}")

    # ---- domain split ------------------------------------------------------
    boundary, score = split_domains(R)
    ecd = [n for n in sorted(R) if n < boundary]
    tmd = [n for n in sorted(R) if n >= boundary]
    print(f"\nECD/TMD boundary computed at residue {boundary} "
          f"(normalised cross-contacts {score:.2f})")
    print(f"  ECD {min(ecd)}-{max(ecd)}  ({len(ecd)} res)")
    print(f"  TMD {min(tmd)}-{max(tmd)}  ({len(tmd)} res)")

    # ---- GLP-1 contacts ----------------------------------------------------
    pc = contacts(P, R, CONTACT_CUTOFF)
    msg_res = [n for n in sorted(P) if n <= MESSAGE_MAX]
    adr_res = [n for n in sorted(P) if n > MESSAGE_MAX]

    msg_hot = sorted({r for p in msg_res for r in pc.get(p, {})})
    adr_hot = sorted({r for p in adr_res for r in pc.get(p, {})})

    print(f"\nGLP-1 message module ({min(msg_res)}-{max(msg_res)}) contacts "
          f"{len(msg_hot)} receptor residues")
    print(f"GLP-1 address helix  ({min(adr_res)}-{max(adr_res)}) contacts "
          f"{len(adr_hot)} receptor residues")

    msg_in_tmd = [r for r in msg_hot if r >= boundary]
    adr_in_ecd = [r for r in adr_hot if r < boundary]
    print(f"  message contacts that fall in the TMD: {len(msg_in_tmd)}/{len(msg_hot)}")
    print(f"  address contacts that fall in the ECD: {len(adr_in_ecd)}/{len(adr_hot)}")

    # ---- cross-check against the literature, do not reconcile --------------
    print("\ncross-check vs brief section 2 (NOT used to define anything):")
    for resnum, label in sorted(LIT_DEEP_POCKET.items()):
        present = resnum in R
        contacted = resnum in msg_hot
        note = ("contacts message" if contacted else
                "contacts address" if resnum in adr_hot else
                "NO contact < %.1f A" % CONTACT_CUTOFF)
        print(f"  {label:>8} {resnum:>4}  {'modelled' if present else 'NOT MODELLED':>12}  {note}")
    print(f"  {'ECL2':>8} {LIT_ECL2:>4}  "
          f"{'modelled' if LIT_ECL2 in R else 'NOT MODELLED':>12}  "
          f"{'contacts GLP-1' if LIT_ECL2 in msg_hot or LIT_ECL2 in adr_hot else 'no contact'}")

    # ---- write targets -----------------------------------------------------
    os.makedirs("target", exist_ok=True)
    n_full = write_pdb(SRC, "R", sorted(R), "target/glp1r_R_full.pdb", "A")
    n_tmd = write_pdb(SRC, "R", tmd, "target/glp1r_TMD.pdb", "A")
    n_ecd = write_pdb(SRC, "R", ecd, "target/glp1r_ECD.pdb", "A")
    print(f"\nwrote target/glp1r_R_full.pdb  {n_full} atoms")
    print(f"wrote target/glp1r_TMD.pdb     {n_tmd} atoms")
    print(f"wrote target/glp1r_ECD.pdb     {n_ecd} atoms")

    out = {
        "source": SRC,
        "source_description": "PDB 6X18 chain R, cryo-EM 2.10 A, active state (Gs + Nb35 bound)",
        "contact_cutoff_A": CONTACT_CUTOFF,
        "domain_boundary_residue": int(boundary),
        "domain_boundary_method": "min normalised cross-domain Ca contacts (<10 A, |i-j|>4)",
        "ecd_range": [int(min(ecd)), int(max(ecd))],
        "tmd_range": [int(min(tmd)), int(max(tmd))],
        "glp1_message_residues": [int(x) for x in msg_res],
        "glp1_address_residues": [int(x) for x in adr_res],
        "targets": {
            "A_orthosteric_vestibule": {
                "pdb": "target/glp1r_TMD.pdb",
                "hotspots": [int(x) for x in msg_in_tmd],
                "rationale": "receptor residues contacted by the GLP-1 message module "
                             "(res 7-15) and lying in the TMD. Class A epitope.",
            },
            "C_ecd_groove": {
                "pdb": "target/glp1r_ECD.pdb",
                "hotspots": [int(x) for x in adr_in_ecd],
                "rationale": "receptor residues contacted by the GLP-1 address helix "
                             "(res 16-36) and lying in the ECD. Class C epitope.",
            },
        },
        "per_glp1_residue_contacts": {
            str(p): {str(r): d for r, d in sorted(pc.get(p, {}).items())}
            for p in sorted(P)
        },
    }
    with open("target/hotspots.json", "w") as f:
        json.dump(out, f, indent=2)
    print("wrote target/hotspots.json")

    print("\nclass A hotspots (TMD):", msg_in_tmd)
    print("class C hotspots (ECD):", adr_in_ecd)


if __name__ == "__main__":
    main()
