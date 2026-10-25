#!/usr/bin/env python3
"""
Build the Boltz-2 shortlist from triage output.

Per-class rules, all from measured calibration rather than imported thresholds:

  class C  AF2-ig is a valid gate (soluble ECD, AF2's home ground).
           Require pae_interaction < AF2_CUT. Rank by pae_interaction, with
           MPNN score breaking ties.

  class A  AF2-ig is NOT a gate: its pae_interaction saturates at ~27.5 for
  class B  every TMD design, which is "no information" (AF2-ig median DockQ
           0.03 on GPCR complexes, brief section 3). Rank on what we can
           actually measure cheaply -- MPNN score, hotspot engagement, contact
           count -- and let Boltz-2 do the discrimination.

Diversity: at most MAX_PER_BACKBONE sequences from any one backbone. Four MPNN
samples of the same scaffold are ~72% identical to each other (measured), so
without this the shortlist collapses onto a handful of backbones and the final
20 cannot satisfy the brief's 70%-identity rule.

Writes <out>/boltz_shortlist.csv (design,sequence) plus a reasons column for
the rejected-candidate log the brief asks for.
"""

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

# AF2-ig gate for class C only. Pre-filter to keep the Boltz load affordable,
# NOT the arbiter -- Boltz-2 is. Calibrated on negative controls, see
# runlog/04-af2-triage-calibration.md.
# Re-calibrated against folded non-binders in the MATCHED representation
# (comp-02 designs at 75/85 aa threaded against the ECD): they score 25.80-27.65.
# Our designs reach 4.35. Cut set well below the null; 145 of 607 pass.
AF2_CUT = 20.0

# Boltz-2 cuts, calibrated against FOLDED NON-BINDERS (comp-02 TNF designs
# scored against GLP1R), which is the null that matters. See
# runlog/05-boltz-calibration.md.
#
#                      POS native   classC best   folded non-binder
#   pae_interaction          4.56         12.55       20.04 - 28.66   <- separates
#   ipTM                    0.949         0.942       0.113 - 0.561   <- separates
#   plddt_binder             91.6          84.1         75.9 - 86.8   <- OVERLAPS
#
# plddt_binder measures foldedness, not binding: our designs sit inside the
# non-binder band. An earlier calibration against SHUFFLED sequences gave the
# opposite answer, because a shuffle tests "does not fold" rather than "does
# not bind". Primary ranking is pae_interaction.
#
# ipTM is used here despite brief section 3 and the comp-02 non-negotiable
# advising against it. Both were reasoning about a HOMOTRIMER target, where
# global ipTM is dominated by target-internal chain pairs at ~0.9. GLP1R is a
# single chain, so ipTM IS the binder-target pair and the separation is wide.
BOLTZ_PAE_CUT = 20.0     # below every folded non-binder
BOLTZ_IPTM_CUT = 0.56    # above every folded non-binder
BOLTZ_PLDDT_MIN = 70.0   # sanity check only: did the binder fold at all

MAX_PER_BACKBONE = 2


def read_csv(path):
    return list(csv.DictReader(open(path))) if os.path.isfile(path) else []


def read_af2(path):
    """AF2-ig writes a Rosetta-style score file, not CSV."""
    out = {}
    if not os.path.isfile(path):
        return out
    hdr = None
    for line in open(path):
        if not line.startswith("SCORE:"):
            continue
        f = line.split()[1:]
        if f[-1] == "description":
            hdr = f
            continue
        if hdr is None:
            continue
        d = dict(zip(hdr, f))
        out[d["description"]] = d
    return out


def cls_of(name):
    for c in ("classA", "classB", "classC"):
        if name.startswith(c):
            return c[-1]
    return "?"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--triage", default="gen/triage")
    p.add_argument("--gen", default="gen/rfdiff")
    p.add_argument("--out", default="gen/triage")
    p.add_argument("--n-c", type=int, default=60)
    p.add_argument("--n-a", type=int, default=40)
    p.add_argument("--n-b", type=int, default=100)
    args = p.parse_args()

    mpnn = read_csv(os.path.join(args.triage, "mpnn/mpnn_manifest.csv"))
    if not mpnn:
        sys.exit(f"no MPNN manifest under {args.triage}")
    af2 = read_af2(os.path.join(args.triage, "af2.sc"))
    denovo = {r["design"]: r for r in read_csv(os.path.join(args.triage,
                                                            "denovo_check.csv"))}

    qc = {}
    for f in glob.glob(os.path.join(args.gen, "interface_qc_*.csv")):
        for r in read_csv(f):
            qc[r["design"]] = r

    rows, rejected = [], []
    for r in mpnn:
        name, parent = r["design"], r["parent_backbone"]
        c = cls_of(name)
        q = qc.get(parent, {})
        dn = denovo.get(name, {})

        if dn and dn.get("denovo_ok", "True") == "False":
            rejected.append((name, "de novo violation: "
                                   f"{dn.get('worst_natural_len')}-residue match to "
                                   f"{dn.get('worst_natural_source')}"))
            continue

        rec = {
            "design": name, "parent_backbone": parent, "class": c,
            "sequence": r["sequence"], "length": r["length"],
            "mpnn_score": float(r["mpnn_score"]) if r.get("mpnn_score") else 9.0,
            "n_contact_res": int(q.get("n_contact_res", 0)),
            "hotspots_engaged": int(q.get("hotspots_engaged", 0)),
            "lipid_frac": float(q.get("lipid_frac", 0) or 0),
        }

        a = af2.get(name)
        if a:
            rec["af2_pae_interaction"] = float(a["pae_interaction"])
            rec["af2_plddt_binder"] = float(a["plddt_binder"])
        else:
            rec["af2_pae_interaction"] = ""
            rec["af2_plddt_binder"] = ""

        if c == "C":
            if not a:
                rejected.append((name, "class C with no AF2 score"))
                continue
            if rec["af2_pae_interaction"] >= AF2_CUT:
                rejected.append((name, f"AF2 pae_interaction "
                                       f"{rec['af2_pae_interaction']:.2f} >= {AF2_CUT}"))
                continue
            rec["rank_key"] = (rec["af2_pae_interaction"], rec["mpnn_score"])
        else:
            # AF2 cannot discriminate here; rank on cheap measurable proxies
            rec["rank_key"] = (rec["mpnn_score"],
                               -rec["hotspots_engaged"],
                               -rec["n_contact_res"])
        rows.append(rec)

    want = {"A": args.n_a, "B": args.n_b, "C": args.n_c}
    shortlist = []
    for c in ("A", "B", "C"):
        pool = sorted([r for r in rows if r["class"] == c],
                      key=lambda r: r["rank_key"])
        per_bb = defaultdict(int)
        taken = 0
        for r in pool:
            if taken >= want[c]:
                break
            if per_bb[r["parent_backbone"]] >= MAX_PER_BACKBONE:
                rejected.append((r["design"], "diversity: backbone already "
                                              f"contributed {MAX_PER_BACKBONE}"))
                continue
            per_bb[r["parent_backbone"]] += 1
            shortlist.append(r)
            taken += 1
        print(f"  class {c}: {len(pool):>5} candidates -> {taken} shortlisted "
              f"from {len(per_bb)} distinct backbones")

    os.makedirs(args.out, exist_ok=True)
    sl = os.path.join(args.out, "boltz_shortlist.csv")
    with open(sl, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "sequence"])
        for r in shortlist:
            w.writerow([r["design"], r["sequence"]])

    full = os.path.join(args.out, "shortlist_full.csv")
    if shortlist:
        keys = [k for k in shortlist[0] if k != "rank_key"]
        with open(full, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(shortlist)

    rej = os.path.join(args.out, "rejected_log.csv")
    with open(rej, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "reason"])
        w.writerows(rejected)

    print(f"\n{len(shortlist)} shortlisted -> {sl}")
    print(f"{len(rejected)} rejected    -> {rej}")
    est = len(shortlist) * 2 * 26.1 / 3600
    print(f"Boltz cost at 2 seeds, 26.1 s/design: {est:.1f} h")


if __name__ == "__main__":
    main()
