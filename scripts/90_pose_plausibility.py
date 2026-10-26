#!/usr/bin/env python3
"""
Is the predicted binder position physically plausible for THIS assay?

Confidence scores and pose reproducibility say nothing about whether a pose is
physically possible in the experiment. Three checks that are specific to this
target and assay format:

1. MEMBRANE PENETRATION. GLP1R sits in a DDM/CHS micelle. A soluble designed
   binder cannot bury itself in the lipid/detergent phase. The membrane slab is
   estimated from the TMD geometry (the hydrophobic belt along the bundle axis),
   and we measure how many binder atoms fall inside it. A binder reaching E387
   at 35 A depth is only credible if it reaches down the aqueous vestibule, not
   through the lipid.

2. C-TERMINUS ACCESSIBILITY. The design is IMMOBILISED on the chip via a
   C-terminal Twin-Strep tag, with receptor flowing over in solution. If the
   C-terminus is buried in the binding interface, immobilising it would
   sterically block the very interaction being measured. The C-terminus must
   point away from the interface.

3. INTERPENETRATION. Severe clashes indicate a physically impossible pose that
   a confidence score will happily rate highly.

Reports per design and flags anything implausible.
"""

import argparse
import csv
import glob
import os

import numpy as np


def load(pdb):
    ch = {}
    for l in open(pdb):
        if l[:4] != "ATOM":
            continue
        if (l[76:78].strip() or l[12:16].strip()[:1]) == "H":
            continue
        ch.setdefault(l[21], {}).setdefault(int(l[22:26]), []).append(
            (l[12:16].strip(),
             np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return ch


def membrane_slab(tmd_ca):
    """
    Estimate the membrane slab from the TM bundle.

    Axis = first principal component of TMD Ca. The hydrophobic belt is taken as
    the central 30 A along that axis, centred on the bundle's midpoint, which is
    the standard thickness for a GPCR-containing bilayer/micelle core.
    """
    X = np.array(list(tmd_ca.values()))
    c = X.mean(0)
    _, _, vt = np.linalg.svd(X - c)
    axis = vt[0] / np.linalg.norm(vt[0])
    proj = (X - c) @ axis
    mid = (proj.min() + proj.max()) / 2
    return c, axis, mid - 15.0, mid + 15.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ranked", default="gen/final_ranked.csv")
    p.add_argument("--struct", default="submission/structures")
    p.add_argument("--out", default="gen/pose_plausibility.csv")
    args = p.parse_args()

    # receptor numbering in the Boltz output is sequential; map to original
    seqmap, i = {}, 0
    for l in open("target/glp1r_R_full.pdb"):
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            i += 1
            seqmap[i] = int(l[22:26])
    inv = {v: k for k, v in seqmap.items()}

    rows = list(csv.DictReader(open(args.ranked)))
    out = []
    for r in rows:
        name = r["design"]
        f = os.path.join(args.struct, f"{name}_complex.pdb")
        if not os.path.isfile(f):
            continue
        ch = load(f)
        binder, target = ch["A"], ch["B"]

        tmd_ca = {n: a[0][1] for n, a in target.items()
                  if 136 <= seqmap.get(n, n) <= 423
                  for an, x in a if an == "CA"}
        tmd_ca = {}
        for n, atoms in target.items():
            if not (136 <= seqmap.get(n, n) <= 423):
                continue
            for an, x in atoms:
                if an == "CA":
                    tmd_ca[n] = x
                    break
        if len(tmd_ca) < 50:
            continue
        c, axis, lo, hi = membrane_slab(tmd_ca)

        batoms = np.array([x for a in binder.values() for _, x in a])
        bproj = (batoms - c) @ axis
        # which side is extracellular: the side the ECD centroid is on
        ecd = [x for n, a in target.items() if seqmap.get(n, n) < 130
               for _, x in a]
        sign = 1.0
        if ecd:
            if ((np.array(ecd) - c) @ axis).mean() < 0:
                sign = -1.0
        in_slab = int(((bproj > lo) & (bproj < hi)).sum())
        frac_slab = in_slab / len(batoms)

        # C-terminus burial: distance from the last residue to the receptor
        cterm = max(binder)
        cxyz = np.array([x for _, x in binder[cterm]])
        tall = np.array([x for a in target.values() for _, x in a])
        d_ct = float(np.linalg.norm(cxyz[:, None, :] - tall[None, :, :],
                                    axis=-1).min())

        # N-terminus, for contrast (the activating element should be here)
        nterm = min(binder)
        nxyz = np.array([x for _, x in binder[nterm]])
        d_nt = float(np.linalg.norm(nxyz[:, None, :] - tall[None, :, :],
                                    axis=-1).min())

        d = np.linalg.norm(batoms[:, None, :] - tall[None, :, :], axis=-1)
        clashes = int((d < 2.2).sum())
        contacts = int((d < 4.5).sum())

        flags = []
        if frac_slab > 0.25:
            flags.append(f"MEMBRANE({100*frac_slab:.0f}%)")
        if d_ct < 5.0:
            flags.append(f"C-TERM-BURIED({d_ct:.1f}A)")
        if clashes > 20:
            flags.append(f"CLASH({clashes})")
        out.append({
            "design": name, "class": r["class"],
            "frac_atoms_in_membrane_slab": round(frac_slab, 3),
            "cterm_min_dist_to_receptor": round(d_ct, 2),
            "nterm_min_dist_to_receptor": round(d_nt, 2),
            "heavy_atom_contacts_4.5A": contacts,
            "clashes_2.2A": clashes,
            "flags": "|".join(flags) or "ok",
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    print(f"{'design':<22}{'cls':>4}{'in-membrane':>13}{'Cterm d':>9}"
          f"{'Nterm d':>9}{'contacts':>10}{'clash':>7}  flags")
    for d in out:
        print(f"{d['design']:<22}{d['class']:>4}"
              f"{100*d['frac_atoms_in_membrane_slab']:>12.0f}%"
              f"{d['cterm_min_dist_to_receptor']:>9.1f}"
              f"{d['nterm_min_dist_to_receptor']:>9.1f}"
              f"{d['heavy_atom_contacts_4.5A']:>10}{d['clashes_2.2A']:>7}"
              f"  {d['flags']}")
    bad = [d for d in out if d["flags"] != "ok"]
    print(f"\n  {len(out)-len(bad)}/{len(out)} plausible; {len(bad)} flagged")
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
