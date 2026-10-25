#!/usr/bin/env python3
"""
Build the final submission package.

Column spec confirmed from the competition page FAQ 5 (see
runlog/06-submission-spec.md), not assumed:

    name            unique identifier
    sequence        the designed protein sequence
    molecule_class  one of protein | nanobody | scfv | fab_kappa | fab_lambda

CSV must be **ordered by rank, top row highest**. Methods, metrics and
structure models are optional but encouraged, so we emit them alongside.

Submitted sequences carry **no Twin-Strep tag** (Harjas, 2026-10-25): the tag is
assay infrastructure appended by the organisers. Every design therefore has to
tolerate a C-terminal fusion, which ours do — the functional element is at the
N-terminus.

Hard gates enforced here, any of which aborts the build:
  * 10 <= length <= 250
  * at most 20 designs
  * single chain, 20 standard amino acids only
  * no contiguous run matching a natural incretin peptide
  * no two designs above 70% identity
  * at least 3 distinct backbone topologies across classes A and B

Writes:
  submission/glp1r_designs.csv     the file to upload
  submission/glp1r_metrics.csv     per-design evidence
  submission/rationale.md          epitope, mechanism, KD evidence, failure mode
  submission/structures/           predicted complexes for the shortlist
"""

import argparse
import csv
import os
import shutil
import sys
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

STD = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 10, 250
MAX_DESIGNS = 20
MAX_IDENTITY = 70.0
MAX_NATURAL_RUN = 5

CLASS_DESC = {
    "A": ("orthosteric vestibule rim (TMD)",
          "Folded miniprotein against the extracellular rim of the GLP-1 "
          "message epitope. Binding only: the activation trigger E387 sits at "
          "35.1 A and is out of reach of a folded scaffold.",
          "thin 6-residue epitope; may not reach high affinity"),
    "B": ("deep orthosteric pocket, two-domain (TMD + ECD)",
          "Extended binder reaching the deep activation-trigger residues "
          "(E387 at 35.1 A, Y241, T391) discovered de novo from the receptor "
          "pocket alone. Intended mechanism: occupy the native agonist pocket "
          "and present an N-terminal element to E387, mimicking the function "
          "of GLP-1's message module without reusing any part of it.",
          "agonism unverified; AF2 cannot score this epitope, Boltz only"),
    "C": ("ECD peptide-binding groove",
          "Miniprotein against the extracellular-domain groove that binds "
          "GLP-1's C-terminal address helix. This is the competition's "
          "recommended epitope.",
          "expected to bind but not activate; published ECD binders are "
          "antagonists"),
}


# Natural incretin/glucagon-family peptides a design must not reproduce.
# Duplicated from 20_denovo_check.py deliberately: this is the last gate before
# upload and it must not depend on a separate script having been run.
NATURAL = {
    "GLP-1(7-36)":   "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
    "GLP-1(1-37)":   "HDEFERHAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG",
    "exendin-4":     "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
    "glucagon":      "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
    "GIP(1-42)":     "YAEGTFISDYSIAMDKIHQQDFVNWLLAQKGKKNDWKHNITQ",
    "oxyntomodulin": "HSQGTFTSDYSKYLDSRRAQDFVQWLMNTKRNRNNIA",
}


def longest_common_substring(a, b):
    best, best_end = 0, 0
    prev = [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best, best_end = cur[j], i
        prev = cur
    return best, a[best_end - best:best_end]


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


def cls_of(name):
    for c in "ABC":
        if name.startswith("class" + c):
            return c
    return "?"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ranked", required=True,
                   help="CSV of final designs, already in rank order, with at "
                        "least design,sequence columns")
    p.add_argument("--out", default="submission")
    p.add_argument("--structures", default="",
                   help="dir of predicted complex PDBs to copy")
    p.add_argument("--force", action="store_true",
                   help="write the package even if gates fail (audit only)")
    args = p.parse_args()

    rows = list(csv.DictReader(open(args.ranked)))
    if not rows:
        sys.exit("ranked file is empty")
    os.makedirs(args.out, exist_ok=True)

    designs = [(r["design"], r["sequence"].strip().upper(), r) for r in rows]
    problems = []

    if len(designs) > MAX_DESIGNS:
        problems.append(f"{len(designs)} designs > cap of {MAX_DESIGNS}")

    for name, seq, _ in designs:
        for src, nat in NATURAL.items():
            L, sub = longest_common_substring(seq, nat)
            if L > MAX_NATURAL_RUN:
                problems.append(f"{name}: DE NOVO VIOLATION, {L}-residue run "
                                f"{sub!r} matches {src}")
        if "HAEG" in seq:
            problems.append(f"{name}: contains HAEG (GLP-1 message motif)")
        if not (MIN_LEN <= len(seq) <= MAX_LEN):
            problems.append(f"{name}: length {len(seq)} outside {MIN_LEN}-{MAX_LEN}")
        bad = set(seq) - STD
        if bad:
            problems.append(f"{name}: non-standard residues {sorted(bad)}")

    for (n1, s1, _), (n2, s2, _) in combinations(designs, 2):
        idp = pct_identity(s1, s2)
        if idp > MAX_IDENTITY:
            problems.append(f"{n1} vs {n2}: {idp:.1f}% identity > {MAX_IDENTITY}")

    # topology diversity across A and B, proxied by distinct parent backbones
    ab_backbones = {r.get("parent_backbone", n) for n, s, r in designs
                    if cls_of(n) in ("A", "B")}
    if len(ab_backbones) < 3:
        problems.append(f"classes A+B span only {len(ab_backbones)} backbones, "
                        f"need >= 3 topologies")

    print(f"{len(designs)} designs")
    if problems:
        print("\nGATE FAILURES:")
        for x in problems:
            print(f"  ! {x}")
        if not args.force:
            sys.exit("\naborting; fix or re-run with --force for an audit build")
    else:
        print("all gates pass")

    # ---- the file that gets uploaded -----------------------------------
    sub = os.path.join(args.out, "glp1r_designs.csv")
    with open(sub, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "sequence", "molecule_class"])
        for name, seq, _ in designs:
            w.writerow([name, seq, "protein"])
    print(f"\nwrote {sub}  ({len(designs)} rows, rank order preserved)")

    # ---- metrics -------------------------------------------------------
    met = os.path.join(args.out, "glp1r_metrics.csv")
    keys = ["design", "class", "length", "mpnn_score", "n_contact_res",
            "hotspots_engaged", "lipid_frac", "af2_pae_interaction",
            "boltz_pae_interaction", "boltz_iptm", "boltz_plddt_binder"]
    with open(met, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for name, seq, r in designs:
            d = dict(r)
            d["design"], d["class"], d["length"] = name, cls_of(name), len(seq)
            w.writerow(d)
    print(f"wrote {met}")

    # ---- per-design rationale -----------------------------------------
    rat = os.path.join(args.out, "rationale.md")
    with open(rat, "w") as f:
        f.write("# Per-design rationale\n\n")
        f.write("Submitted sequences carry no Twin-Strep tag; the C-terminal "
                "capture tag is assay infrastructure added by the organisers. "
                "Every design places its functional element at the N-terminus "
                "and tolerates a C-terminal fusion.\n\n")
        for i, (name, seq, r) in enumerate(designs, 1):
            c = cls_of(name)
            ep, mech, fail = CLASS_DESC.get(c, ("?", "?", "?"))
            f.write(f"## {i}. {name}\n\n")
            f.write(f"- **length**: {len(seq)} aa\n")
            f.write(f"- **epitope targeted**: {ep}\n")
            f.write(f"- **intended mechanism**: {mech}\n")
            bp = r.get("boltz_pae_interaction", "")
            bi = r.get("boltz_iptm", "")
            f.write(f"- **predicted-KD evidence**: Boltz-2 "
                    f"pae_interaction {bp or 'n/a'}, ipTM {bi or 'n/a'}. "
                    f"No KD is predicted; these are interface-confidence "
                    f"proxies calibrated against folded non-binders, not "
                    f"affinity estimates.\n")
            f.write(f"- **main failure mode**: {fail}\n\n")
    print(f"wrote {rat}")

    if args.structures and os.path.isdir(args.structures):
        sd = os.path.join(args.out, "structures")
        os.makedirs(sd, exist_ok=True)
        n = 0
        for name, _, _ in designs:
            for cand in (f"{name}_model_0.pdb", f"{name}_af2pred.pdb"):
                src = os.path.join(args.structures, cand)
                if os.path.isfile(src):
                    shutil.copy(src, sd)
                    n += 1
                    break
        print(f"copied {n} predicted structures -> {sd}")

    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
