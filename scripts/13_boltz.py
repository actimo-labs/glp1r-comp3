#!/usr/bin/env python3
"""
Stage 4 (orthogonal validation): Boltz-2 on binder + GLP1R.

Scoring follows brief section 3 Stage 4: rank on INTER-CHAIN PAE and BINDER
pLDDT, not ipTM. On the 113-complex GPCR benchmark quoted there, Boltz-2's
confidence metrics systematically overestimate correctness for misplaced
peptides, and medium vs high quality were indistinguishable by PAE alone -- so
these numbers are a filter, not a ranking you can trust on its own. ipTM is
recorded but never used to rank.

Target templating is not optional here. Two reasons:
  * comp-02 found that without a template Boltz folds the chains but does not
    dock them, so the designed epitope is simply absent from the prediction.
  * we are designing against the ACTIVE state specifically. Left free, Boltz
    may return some other receptor conformation and the binder pose would then
    be scored against the wrong geometry.
So the receptor is pinned to target/glp1r_R_full.cif (6X18 chain R).

Run under the boltz env:
  ~/miniforge3/envs/boltz/bin/python scripts/13_boltz.py --manifest ... --out ...
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from timeit import default_timer as timer

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TARGET_SEQ = (ROOT / "target/glp1r_R_full.seq").read_text().strip()
TEMPLATE = ROOT / "target/glp1r_R_full.cif"
BINDER_ID, TARGET_ID = "A", "B"


def write_yaml(path, binder_seq):
    lines = ["version: 1", "sequences:",
             "  - protein:", f"      id: {BINDER_ID}",
             f"      sequence: {binder_seq}", "      msa: empty",
             "  - protein:", f"      id: {TARGET_ID}",
             f"      sequence: {TARGET_SEQ}", "      msa: empty",
             "templates:", f"  - cif: {TEMPLATE}",
             f"    chain_id: [{TARGET_ID}]", "    template_id: [A]"]
    path.write_text("\n".join(lines) + "\n")


def parse_prediction(pred_root, name):
    """Return inter-chain PAE and per-chain pLDDT for one design."""
    conf_hits = sorted(Path(pred_root).rglob(f"confidence_{name}_model_0.json"))
    pae_hits = sorted(Path(pred_root).rglob(f"pae_{name}_model_0.npz"))
    pdb_hits = (sorted(Path(pred_root).rglob(f"{name}_model_0.pdb"))
                or sorted(Path(pred_root).rglob(f"{name}_model_0.cif")))
    if not conf_hits or not pdb_hits:
        return None

    conf = json.loads(conf_hits[0].read_text())

    # chain membership + per-residue pLDDT from the output structure B-factors
    chains, plddt = [], []
    seen = set()
    for l in pdb_hits[0].read_text().splitlines():
        if l[:4] != "ATOM":
            continue
        if l[12:16].strip() != "CA":
            continue
        key = (l[21], l[22:27])
        if key in seen:
            continue
        seen.add(key)
        chains.append(l[21])
        try:
            plddt.append(float(l[60:66]))
        except ValueError:
            plddt.append(float("nan"))
    chains = np.array(chains)
    plddt = np.array(plddt, dtype=float)
    if plddt.max() <= 1.01:       # boltz writes 0-1 in cif, 0-100 in pdb
        plddt = plddt * 100.0
    is_b = chains == BINDER_ID
    is_t = ~is_b

    row = {
        "design": name,
        "plddt_binder": float(np.nanmean(plddt[is_b])) if is_b.any() else float("nan"),
        "plddt_target": float(np.nanmean(plddt[is_t])) if is_t.any() else float("nan"),
        "plddt_complex": float(np.nanmean(plddt)),
        "binderlen": int(is_b.sum()),
        # recorded, deliberately NOT used to rank -- see module docstring
        "iptm_recorded": conf.get("iptm"),
        "ptm_recorded": conf.get("ptm"),
    }

    if pae_hits:
        pae = np.load(pae_hits[0])
        pae = pae[pae.files[0]]
        if pae.ndim == 3:
            pae = pae[0]
        if pae.shape[0] == len(chains):
            bt = pae[np.ix_(is_b, is_t)]
            tb = pae[np.ix_(is_t, is_b)]
            row["pae_interaction"] = float((bt.mean() + tb.mean()) / 2)
            row["pae_binder"] = float(pae[np.ix_(is_b, is_b)].mean())
            row["pae_interaction_min"] = float(min(bt.mean(), tb.mean()))
        else:
            row["pae_interaction"] = float("nan")
    else:
        row["pae_interaction"] = float("nan")
    return row


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, help="CSV with design,sequence")
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--recycling", type=int, default=3)
    p.add_argument("--diffusion-samples", type=int, default=1)
    p.add_argument("--limit", type=int, default=0)
    args = p.parse_args()

    if not TEMPLATE.is_file():
        sys.exit(f"missing template {TEMPLATE}; run make_template_cif.py")

    rows = list(csv.DictReader(open(args.manifest)))
    if args.limit:
        rows = rows[:args.limit]
    out = Path(args.out)
    yaml_dir = out / "yaml"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    for r in rows:
        write_yaml(yaml_dir / f"{r['design']}.yaml", r["sequence"].strip())
    print(f"[boltz] {len(rows)} yaml -> {yaml_dir}")

    t0 = timer()
    cmd = ["boltz", "predict", str(yaml_dir),
           "--out_dir", str(out),
           "--model", "boltz2",
           "--output_format", "pdb",
           "--write_full_pae",
           "--recycling_steps", str(args.recycling),
           "--diffusion_samples", str(args.diffusion_samples),
           "--seed", str(args.seed),
           "--accelerator", "gpu", "--devices", "1",
           "--override"]
    print("[boltz]", " ".join(cmd))
    rc = subprocess.run(cmd).returncode
    wall = timer() - t0
    print(f"[boltz] rc={rc}  {wall:.1f}s  "
          f"({wall/max(1,len(rows)):.1f}s/design)")

    parsed, missing = [], []
    for r in rows:
        got = parse_prediction(out, r["design"])
        if got is None:
            missing.append(r["design"])
            continue
        got["sequence"] = r["sequence"].strip()
        parsed.append(got)

    if parsed:
        scorefile = out / "boltz_scores.csv"
        keys = sorted({k for d in parsed for k in d})
        with open(scorefile, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(parsed)
        parsed.sort(key=lambda d: (d.get("pae_interaction") or 99))
        print(f"\n{'design':<24}{'pae_int':>9}{'plddt_bind':>12}{'iptm(rec)':>11}")
        for d in parsed[:20]:
            print(f"{d['design']:<24}{d.get('pae_interaction', float('nan')):>9.2f}"
                  f"{d.get('plddt_binder', float('nan')):>12.1f}"
                  f"{(d.get('iptm_recorded') or float('nan')):>11.3f}")
        print(f"\nwrote {scorefile}")
    if missing:
        print(f"[boltz] {len(missing)} with no parsable output: {missing[:10]}")
    return 0 if parsed else 1


if __name__ == "__main__":
    sys.exit(main())
