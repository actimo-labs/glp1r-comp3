#!/usr/bin/env python3
"""
Final ranking and the 20-design cut.

Gates, in order. A design must clear all of them:

  1. self-consistency   binder RMSD < 5 A between the Boltz prediction and the
                        RFdiffusion design, after superposing on the receptor.
                        Boltz re-predicts the binder pose, so a good interface
                        score on a pose tens of angstroms from the design is not
                        a validated design -- it is a different complex. This is
                        the gate that eliminated class B entirely (median 64.8 A).
  2. Boltz interface    pae_interaction < 20.0 and ipTM > 0.56, both calibrated
                        against folded non-binders (20.04-28.66 / 0.113-0.561).
  3. epitope            interface not predominantly on the lipid-facing TM
                        surface (RFdiffusion has no membrane; this is its
                        favourite failure mode on a GPCR).
  4. de novo            no contiguous run > 5 residues matching a natural
                        incretin peptide.

Ranking among survivors: pae_interaction ascending, ipTM descending as a
tie-break. plddt_binder is NOT used -- it was shown three times over to track
foldedness rather than binding (folded non-binders score 83-94, the same band
as our best designs).

Diversity: at most 2 designs per parent backbone, and no pair above 70%
identity, enforced greedily down the ranking.
"""

import argparse
import csv
import glob
import json
import os
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np

RMSD_CUT = 5.0
PAE_CUT = 20.0
IPTM_CUT = 0.56
LIPID_CUT = 0.35
MAX_PER_BACKBONE = 2
MAX_IDENTITY = 70.0

TRIGGER = {387: "E387", 234: "Q234", 233: "M233", 391: "T391", 241: "Y241"}


def pct_identity(a, b, min_overlap_frac=0.8):
    if not a or not b:
        return 0.0
    shorter = min(len(a), len(b))
    need = max(1, int(shorter * min_overlap_frac))
    best = 0
    for off in range(-(len(b) - 1), len(a)):
        m = n = 0
        for i in range(len(a)):
            j = i - off
            if 0 <= j < len(b):
                n += 1
                m += a[i] == b[j]
        if n >= need:
            best = max(best, m / shorter)
    return best * 100


def cls_of(n):
    return "A" if n.startswith("classA") else "B" if n.startswith("classB") else "C"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sc", default="gen/boltz/selfconsistency.csv")
    p.add_argument("--manifest", default="gen/triage/mpnn_all.csv")
    p.add_argument("--denovo", default="gen/triage/denovo_all.csv")
    p.add_argument("--boltzdir", default="gen/boltz")
    p.add_argument("--out", default="gen/final_ranked.csv")
    p.add_argument("--n", type=int, default=20)
    # Per-class quotas rather than a single global ranking. Ranking purely on
    # pae_interaction put 19 of 20 slots on class A, which concentrates the
    # whole submission on one epitope. Two reasons not to:
    #   - the competition states that binding at the peptide groove (the ECD,
    #     i.e. class C) counts in a design's favour
    #   - if class A's predicted deep insertion is a Boltz artefact, a
    #     single-epitope submission fails entirely. Splitting hedges that.
    p.add_argument("--n-a", type=int, default=10)
    p.add_argument("--n-c", type=int, default=10)
    args = p.parse_args()

    ec = json.load(open("target/epitope_classes.json"))
    lipid = set(ec["tmd_lipid_facing_all"])
    seqmap, i = {}, 0
    for l in open("target/glp1r_R_full.pdb"):
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            i += 1
            seqmap[i] = int(l[22:26])

    seqs = {r["design"]: r["sequence"] for r in csv.DictReader(open(args.manifest))}
    parents = {r["design"]: r["parent_backbone"]
               for r in csv.DictReader(open(args.manifest))}
    mscore = {r["design"]: r.get("mpnn_score", "")
              for r in csv.DictReader(open(args.manifest))}
    dn = {r["design"]: r["denovo_ok"] for r in csv.DictReader(open(args.denovo))} \
        if os.path.isfile(args.denovo) else {}

    rows = list(csv.DictReader(open(args.sc)))
    kept, rejected = [], []

    for r in rows:
        n = r["design"]
        rmsd = float(r["binder_rmsd_vs_design"])
        pae = float(r["pae_interaction"])
        iptm = float(r.get("iptm_recorded") or 0)

        if rmsd >= RMSD_CUT:
            rejected.append((n, f"not self-consistent: binder RMSD {rmsd:.1f} A"))
            continue
        if pae >= PAE_CUT:
            rejected.append((n, f"Boltz pae_interaction {pae:.2f} >= {PAE_CUT}"))
            continue
        if iptm <= IPTM_CUT:
            rejected.append((n, f"Boltz ipTM {iptm:.3f} <= {IPTM_CUT}"))
            continue
        if dn.get(n) == "False":
            rejected.append((n, "de novo violation"))
            continue

        hits = glob.glob(os.path.join(args.boltzdir, "**", f"{n}_model_0.pdb"),
                         recursive=True)
        contacts, lf, trig = set(), 0.0, []
        if hits:
            ch = {}
            for l in open(hits[0]):
                if l[:4] == "ATOM":
                    ch.setdefault(l[21], {}).setdefault(int(l[22:26]), []).append(
                        [float(l[30:38]), float(l[38:46]), float(l[46:54])])
            if "A" in ch and "B" in ch:
                b = [np.array(v) for v in ch["A"].values()]
                for tn, tv in ch["B"].items():
                    tv = np.array(tv)
                    if min(np.linalg.norm(x[:, None, :] - tv[None, :, :],
                                          axis=-1).min() for x in b) < 5.0:
                        contacts.add(seqmap.get(tn, tn))
                lf = sum(1 for c in contacts if c in lipid) / max(1, len(contacts))
                trig = [v for k, v in TRIGGER.items() if k in contacts]
        if lf > LIPID_CUT:
            rejected.append((n, f"interface {100*lf:.0f}% lipid-facing"))
            continue

        kept.append({
            "design": n, "class": cls_of(n), "sequence": seqs.get(n, ""),
            "length": len(seqs.get(n, "")),
            "parent_backbone": parents.get(n, ""),
            "boltz_pae_interaction": round(pae, 2),
            "boltz_iptm": round(iptm, 3),
            "boltz_plddt_binder": round(float(r.get("plddt_binder") or 0), 1),
            "binder_rmsd_vs_design": rmsd,
            "n_contacts_pred": len(contacts),
            "lipid_frac": round(lf, 3),
            "trigger_contacts": "|".join(trig),
            "n_trigger": len(trig),
            "mpnn_score": mscore.get(n, ""),
        })

    kept.sort(key=lambda d: (d["boltz_pae_interaction"], -d["boltz_iptm"]))

    quota = {"A": args.n_a, "C": args.n_c, "B": 0}
    taken = defaultdict(int)
    final, per_bb = [], defaultdict(int)
    for d in kept:
        if len(final) >= args.n:
            break
        if taken[d["class"]] >= quota.get(d["class"], 0):
            rejected.append((d["design"], f"class {d['class']} quota filled"))
            continue
        if per_bb[d["parent_backbone"]] >= MAX_PER_BACKBONE:
            rejected.append((d["design"], "diversity: backbone quota"))
            continue
        clash = [f for f in final
                 if pct_identity(f["sequence"], d["sequence"]) > MAX_IDENTITY]
        if clash:
            rejected.append((d["design"],
                             f"identity > {MAX_IDENTITY}% vs {clash[0]['design']}"))
            continue
        per_bb[d["parent_backbone"]] += 1
        taken[d["class"]] += 1
        final.append(d)

    if not final:
        sys.exit("nothing survived the gates")

    keys = list(final[0].keys())
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(final)
    rej = args.out.replace(".csv", "_rejected.csv")
    with open(rej, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "reason"])
        w.writerows(rejected)

    byc = defaultdict(int)
    for d in final:
        byc[d["class"]] += 1
    print(f"{len(kept)} designs cleared all gates; top {len(final)} selected")
    print(f"  composition: " + "  ".join(f"class {c}: {byc[c]}" for c in sorted(byc)))
    print(f"  backbones:   {len(per_bb)} distinct")
    print(f"  with E387 contact: {sum(1 for d in final if 'E387' in d['trigger_contacts'])}")
    print(f"\n{'rank':<5}{'design':<22}{'cls':>4}{'pae':>7}{'ipTM':>7}{'rmsd':>6}{'trig':>6}  triggers")
    for i, d in enumerate(final, 1):
        print(f"{i:<5}{d['design']:<22}{d['class']:>4}{d['boltz_pae_interaction']:>7.2f}"
              f"{d['boltz_iptm']:>7.3f}{d['binder_rmsd_vs_design']:>6.1f}"
              f"{d['n_trigger']:>6}  {d['trigger_contacts']}")
    print(f"\nwrote {args.out}\nwrote {rej}")


if __name__ == "__main__":
    main()
