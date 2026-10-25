#!/usr/bin/env python3
"""
State selectivity: does each design's binding site exist in the INACTIVE receptor?

Why not just compute activation order parameters on our predicted complexes:
Boltz pins the receptor to the active-state template (6X18 chain R), so every
prediction shows active geometry no matter what the binder does. Measuring TM6
position on those is circular and proves nothing. The brief makes the same point
about stripping Gs in MD -- the receptor's state is imposed, not discovered.

The non-circular question is a steric one, and it maps onto conformational
selection: an agonist preferentially binds a pre-existing active-like
conformation. So if a design's epitope is only open in the active state, the
design can only engage the active state, which is the structural signature we
want. If it fits both states equally, it has no reason to shift the equilibrium
and is more likely neutral or an antagonist.

Method:
  1. superpose the inactive receptor onto the active one using TMD-core Ca
     atoms common to both (ECD and loops excluded -- the brief notes global RMSD
     is dominated by ECD motion and is the wrong order parameter)
  2. place each design's binder, which is already in the active frame, against
     the superposed inactive receptor
  3. count heavy-atom clashes (< 2.5 A) and lost contacts

  clash_inactive high + clash_active low  -> state-selective, agonist candidate
  clash_inactive ~ clash_active           -> state-agnostic, likely not an agonist

This is a STRUCTURAL PROXY. It is not a prediction of cAMP activation and must
not be reported as one.
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np

CLASH_A = 2.5
CONTACT_A = 5.0


def ca_by_res(pdb, chain=None, lo=None, hi=None):
    out = {}
    for l in open(pdb):
        if l[:4] != "ATOM" or l[12:16].strip() != "CA":
            continue
        if chain and l[21] != chain:
            continue
        if l[16] not in (" ", "A"):
            continue
        n = int(l[22:26])
        if lo is not None and not (lo <= n <= hi):
            continue
        out[n] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


def heavy(pdb, chain=None, lo=None, hi=None):
    out = []
    for l in open(pdb):
        if l[:4] != "ATOM":
            continue
        if chain and l[21] != chain:
            continue
        if l[16] not in (" ", "A"):
            continue
        if (l[76:78].strip() or l[12:16].strip()[:1]) == "H":
            continue
        n = int(l[22:26])
        if lo is not None and not (lo <= n <= hi):
            continue
        out.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(out)


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def pick_receptor_chain(pdb):
    """The GLP1R chain is the longest polypeptide with H363 present."""
    best, bestn = None, 0
    chains = {}
    for l in open(pdb):
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            chains.setdefault(l[21], set()).add(int(l[22:26]))
    for c, res in chains.items():
        if 363 in res and len(res) > bestn:
            best, bestn = c, len(res)
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inactive", required=True)
    p.add_argument("--ranked", default="gen/final_ranked.csv")
    p.add_argument("--boltzdir", default="gen/boltz")
    p.add_argument("--tmd-lo", type=int, default=136)
    p.add_argument("--tmd-hi", type=int, default=423)
    p.add_argument("--out", default="gen/state_selectivity.csv")
    args = p.parse_args()

    ichain = pick_receptor_chain(args.inactive)
    if ichain is None:
        sys.exit(f"{args.inactive}: could not identify the GLP1R chain")
    print(f"inactive {os.path.basename(args.inactive)}: receptor chain {ichain}")

    act_ca = ca_by_res("target/glp1r_R_full.pdb", "A", args.tmd_lo, args.tmd_hi)
    ina_ca = ca_by_res(args.inactive, ichain, args.tmd_lo, args.tmd_hi)
    shared = sorted(set(act_ca) & set(ina_ca))
    print(f"  TMD-core Ca shared with the active structure: {len(shared)}")
    if len(shared) < 50:
        sys.exit("too few shared TMD residues to superpose reliably")

    P = np.array([ina_ca[n] for n in shared])
    Q = np.array([act_ca[n] for n in shared])
    R, t = kabsch(P, Q)
    rms = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
    print(f"  TMD-core superposition RMSD: {rms:.2f} A")

    # sanity: TM6 (K346 Ca) should sit far apart between the two states
    if 346 in act_ca and 346 in ina_ca:
        d = np.linalg.norm((R @ ina_ca[346] + t) - act_ca[346])
        print(f"  TM6 reference K346 Ca displacement active vs inactive: {d:.1f} A")

    ina_heavy = heavy(args.inactive, ichain)
    ina_heavy = (R @ ina_heavy.T).T + t
    act_heavy = heavy("target/glp1r_R_full.pdb", "A")

    rows = list(csv.DictReader(open(args.ranked)))
    out = []
    for r in rows:
        n = r["design"]
        hits = glob.glob(os.path.join(args.boltzdir, "**", f"{n}_model_0.pdb"),
                         recursive=True)
        if not hits:
            continue
        binder = heavy(hits[0], "A")
        if binder.size == 0:
            continue

        def counts(target):
            d = np.linalg.norm(binder[:, None, :] - target[None, :, :], axis=-1)
            return int((d < CLASH_A).sum()), int((d < CONTACT_A).sum())

        ci, cti = counts(ina_heavy)
        ca_, cta = counts(act_heavy)
        sel = (ci - ca_) / max(1, binder.shape[0]) * 100
        rec = dict(r)
        rec.update({
            "clash_active": ca_, "clash_inactive": ci,
            "contact_active": cta, "contact_inactive": cti,
            "clash_per100atom_active": round(100 * ca_ / binder.shape[0], 1),
            "clash_per100atom_inactive": round(100 * ci / binder.shape[0], 1),
            "state_selectivity": round(sel, 1),
        })
        out.append(rec)

    if not out:
        sys.exit("no designs scored")
    out.sort(key=lambda d: -d["state_selectivity"])
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    print(f"\n{'design':<22}{'cls':>4}{'clash_act':>11}{'clash_ina':>11}"
          f"{'selectivity':>13}  interpretation")
    for d in out:
        s = d["state_selectivity"]
        tag = ("active-selective" if s > 5 else
               "state-agnostic" if s > -5 else "inactive-selective")
        print(f"{d['design']:<22}{d['class']:>4}"
              f"{d['clash_per100atom_active']:>11.1f}{d['clash_per100atom_inactive']:>11.1f}"
              f"{s:>13.1f}  {tag}")
    print(f"\nwrote {args.out}")
    print("\nSTRUCTURAL PROXY ONLY -- not a prediction of cAMP activation.")


if __name__ == "__main__":
    main()
