#!/usr/bin/env python3
"""
Interface QC on RFdiffusion backbones, before spending any sequence-design time.

RFdiffusion emits the binder as chain A (residues 1..L) and renumbers the target
into chain B starting at L+1, so target numbering must be mapped back to the
original receptor numbering before hotspots mean anything.

Rejects, in order of how badly they waste downstream compute:

  LIPID    the interface sits on the lipid-facing TM surface. RFdiffusion has no
           membrane, so this is its favourite failure on a GPCR: a big greasy
           convex patch. These cannot bind in a detergent micelle.
  NOCONTACT too few binder-target contacts to be a real interface.
  OFFTARGET contacts exist but miss the requested hotspots entirely.

Prints a table and writes <dir>/interface_qc.csv.
"""

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np

MIN_CONTACT_RES = 8     # target residues within 5 A of the binder
CONTACT_A = 5.0
MAX_LIPID_FRAC = 0.35   # fraction of the interface allowed on lipid-facing surface


def load_chains(pdb):
    ch = {}
    for l in open(pdb):
        if l[:4] != "ATOM":
            continue
        ch.setdefault(l[21], {}).setdefault(int(l[22:26]), []).append(
            [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return {c: {n: np.array(v) for n, v in r.items()} for c, r in ch.items()}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", required=True)
    p.add_argument("--tag", default="", help="only files starting with this")
    p.add_argument("--hotspots", required=True,
                   help="comma-separated, original receptor numbering")
    p.add_argument("--target-start", type=int, required=True,
                   help="first ORIGINAL residue number of the target as passed "
                        "to RFdiffusion (e.g. 29 for the ECD, 136 for the TMD)")
    args = p.parse_args()

    hot = [int(x) for x in args.hotspots.split(",") if x]
    ec = json.load(open("target/epitope_classes.json"))
    lipid = set(ec["tmd_lipid_facing_all"])

    pdbs = sorted(glob.glob(os.path.join(args.dir, f"{args.tag}*.pdb")))
    if not pdbs:
        sys.exit(f"no PDBs matching {args.tag}* in {args.dir}")

    rows = []
    for f in pdbs:
        ch = load_chains(f)
        if "A" not in ch or "B" not in ch:
            continue
        binder, target = ch["A"], ch["B"]
        off = min(target) - args.target_start     # output -> original

        bl = list(binder.values())
        contacts = {}
        for tn, tv in target.items():
            d = min(np.linalg.norm(bv[:, None, :] - tv[None, :, :], axis=-1).min()
                    for bv in bl)
            if d < CONTACT_A:
                contacts[tn - off] = round(float(d), 2)

        engaged = sorted(h for h in hot if h in contacts)
        n_lip = sum(1 for r in contacts if r in lipid)
        lip_frac = n_lip / max(1, len(contacts))

        if len(contacts) < MIN_CONTACT_RES:
            verdict = "NOCONTACT"
        elif lip_frac > MAX_LIPID_FRAC:
            verdict = "LIPID"
        elif not engaged:
            verdict = "OFFTARGET"
        else:
            verdict = "ok"

        rows.append({
            "design": os.path.basename(f)[:-4],
            "binder_len": len(binder),
            "n_contact_res": len(contacts),
            "hotspots_engaged": len(engaged),
            "hotspots_list": "|".join(map(str, engaged)),
            "lipid_frac": round(lip_frac, 2),
            "verdict": verdict,
        })

    rows.sort(key=lambda r: (r["verdict"] != "ok", -r["hotspots_engaged"],
                             -r["n_contact_res"]))
    print(f"{'design':<16}{'len':>5}{'contacts':>10}{'hot':>5}{'lipid':>7}  verdict")
    for r in rows:
        print(f"{r['design']:<16}{r['binder_len']:>5}{r['n_contact_res']:>10}"
              f"{r['hotspots_engaged']:>5}{r['lipid_frac']:>7.2f}  {r['verdict']}"
              + (f"  [{r['hotspots_list']}]" if r["hotspots_list"] else ""))

    out = os.path.join(args.dir, f"interface_qc{('_'+args.tag) if args.tag else ''}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_ok = sum(1 for r in rows if r["verdict"] == "ok")
    print(f"\n{n_ok}/{len(rows)} pass  ->  {out}")
    from collections import Counter
    for v, c in Counter(r["verdict"] for r in rows).most_common():
        print(f"  {v:<10} {c}")


if __name__ == "__main__":
    main()
