#!/usr/bin/env python3
"""
Build two-sided scoring controls.

comp-02's calibration (scripts/04b_make_controls.py) used composition-matched
shuffles as negatives and said plainly what that buys: it calibrates the null,
so the cutoff means "clearly better than noise", not "comparable to a real
binder". It had no known-binder positive available.

This target does. Native GLP-1 is a validated sub-nanomolar binder AND a full
agonist, and its bound structure IS the design target (6X18 chain P). So both
ends can be anchored:

  POSITIVES  native GLP-1(7-36) and exendin-4 on the native backbone.
             These are the ceiling: whatever score a real, validated,
             experimentally-confirmed agonist earns in this setup is the best
             any design can be expected to do. They are controls only and are
             never submitted -- they would fail the de novo rule outright,
             which scripts/20_denovo_check.py confirms.

  NEGATIVES  composition- and length-matched shuffles of the same sequences,
             on the same backbone. Identical amino-acid content, destroyed
             sequence-structure relationship. Whatever these still earn is the
             floor.

Why this matters concretely: the native GLP-1 positive scored plddt_binder 84.0
through AF2-ig here. Brief section 3 proposes Baker's class-B thresholds of
plddt_binder > 85-90. Applied literally that threshold REJECTS the native
agonist, so it cannot be right for this target and metric. Thresholds get set
from these controls instead.

Writes <out>/controls.csv (design,sequence) ready for 13_boltz.py, and
<out>/threaded/*.pdb ready for the AF2-ig driver.
"""

import argparse
import csv
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ONE2THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
BACKBONE = ("N", "CA", "C", "O")

POSITIVES = {
    # GLP-1(7-36)NH2 -- exactly the peptide in 6X18 chain P
    "POS_native_glp1": "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
    # exendin-4 N-terminal 30 -- same length, also a validated full agonist
    "POS_exendin4_30": "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGG",
}


def first_chain(pdb):
    """Binder chain = whichever chain appears first. The convention everywhere
    in this pipeline is binder-first, but the native complex keeps 6X18's own
    chain IDs (P/R) rather than A/B, so this must not be hardcoded."""
    for l in open(pdb):
        if l[:4] == "ATOM":
            return l[21]
    raise ValueError(f"{pdb}: no ATOM records")


def thread(backbone_pdb, seq, out_pdb, binder_chain=None):
    """Same threading rule as 11_mpnn_thread.py: backbone atoms, renamed."""
    if binder_chain is None:
        binder_chain = first_chain(backbone_pdb)
    lines, binder_res = [], []
    for l in open(backbone_pdb):
        if l[:4] != "ATOM":
            continue
        lines.append(l)
        if l[21] == binder_chain:
            n = int(l[22:26])
            if not binder_res or binder_res[-1] != n:
                binder_res.append(n)
    if len(binder_res) != len(seq):
        raise ValueError(f"backbone has {len(binder_res)} binder residues, "
                         f"sequence has {len(seq)}")
    aamap = {n: ONE2THREE[a] for n, a in zip(binder_res, seq)}
    serial = 0
    with open(out_pdb, "w") as f:
        for l in lines:
            if l[21] == binder_chain:
                if l[12:16].strip() not in BACKBONE:
                    continue
                serial += 1
                f.write(l[:6] + f"{serial:5d}" + l[11:17]
                        + aamap[int(l[22:26])] + l[20:])
            else:
                serial += 1
                f.write(l[:6] + f"{serial:5d}" + l[11:])
        f.write("TER\nEND\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--backbone", default="pilot/af2_smoke_in/native_glp1_glp1r.pdb",
                   help="native complex, binder chain first")
    p.add_argument("--out", default="controls")
    p.add_argument("--n-shuffles", type=int, default=4)
    p.add_argument("--seed", type=int, default=17)
    args = p.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(os.path.join(args.out, "threaded"), exist_ok=True)

    rows = []
    for name, seq in POSITIVES.items():
        rows.append((name, seq, "positive"))
        for i in range(args.n_shuffles):
            sl = list(seq)
            rng.shuffle(sl)
            shuf = "".join(sl)
            rows.append((f"NEG_shuf_{name.split('_',1)[1]}_{i}", shuf, "negative"))

    n_ok = 0
    for name, seq, kind in rows:
        try:
            thread(args.backbone, seq,
                   os.path.join(args.out, "threaded", f"{name}.pdb"))
            n_ok += 1
        except Exception as e:
            print(f"  {name}: FAILED {e}")

    man = os.path.join(args.out, "controls.csv")
    with open(man, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "sequence", "control_type"])
        w.writerows(rows)

    print(f"{len(rows)} controls "
          f"({sum(1 for r in rows if r[2]=='positive')} positive, "
          f"{sum(1 for r in rows if r[2]=='negative')} negative), "
          f"{n_ok} threaded")
    print(f"wrote {man}")
    print("\nNOTE: positives are native peptides. They are calibration only and")
    print("must never enter the submission -- 20_denovo_check.py fails them by")
    print("design, which is the intended behaviour.")


if __name__ == "__main__":
    main()
