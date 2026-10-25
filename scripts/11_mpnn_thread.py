#!/usr/bin/env python3
"""
Stage: ProteinMPNN sequence design on RFdiffusion backbones, then thread the
designed sequences back onto those backbones as AF2-initial-guess inputs.

Replaces dl_binder_design's `mpnn_fr` stage, which is MPNN + Rosetta FastRelax
and needs PyRosetta. We have no licence, so there is no FastRelax here. That is
a real omission and is recorded in the methods "not computed / not run" list:
upstream relaxes the backbone between MPNN and AF2, and we do not.

Input  : directory of RFdiffusion outputs. Binder MUST be chain A (RFdiffusion
         emits it first, so this is automatic), target chain B.
Output : <out>/seqs/*.fa          raw MPNN output
         <out>/threaded/*.pdb     AF2-ig-ready, binder chain A with the designed
                                  sequence on the RFdiffusion backbone
         <out>/mpnn_manifest.csv  design -> sequence, score, parent backbone

Threading note: RFdiffusion emits the binder as poly-GLY with N/CA/C/O only and
no CB. We keep those four backbone atoms and rename the residue. AF2's atom
parser maps by name into the 37-atom representation and masks everything it does
not find, which is exactly what "initial guess" wants for the chain being
predicted. The target keeps its full side chains, because it is templated.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import glob
from timeit import default_timer as timer

MPNN = os.path.expanduser("~/tnf-track3/ProteinMPNN")
MPNN_PY = os.path.expanduser("~/miniforge3/envs/mpnn/bin/python")

ONE2THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
BACKBONE = ("N", "CA", "C", "O")


def sh(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.stderr.write(f"FAILED: {' '.join(cmd)}\n{r.stdout}\n{r.stderr}\n")
        sys.exit(r.returncode)
    return r


def run_mpnn(pdbdir, out, nseq, temp, seed, omit):
    os.makedirs(out, exist_ok=True)
    jsonl = os.path.join(out, "parsed.jsonl")
    chains = os.path.join(out, "chains.jsonl")

    sh([MPNN_PY, f"{MPNN}/helper_scripts/parse_multiple_chains.py",
        "--input_path", pdbdir, "--output_path", jsonl])
    # chain A designable, everything else fixed
    sh([MPNN_PY, f"{MPNN}/helper_scripts/assign_fixed_chains.py",
        "--input_path", jsonl, "--output_path", chains, "--chain_list", "A"])

    t0 = timer()
    sh([MPNN_PY, f"{MPNN}/protein_mpnn_run.py",
        "--jsonl_path", jsonl,
        "--chain_id_jsonl", chains,
        "--out_folder", out,
        "--num_seq_per_target", str(nseq),
        "--sampling_temp", str(temp),
        "--seed", str(seed),
        "--batch_size", "1",
        "--omit_AAs", omit,
        "--use_soluble_model"])
    return timer() - t0


def parse_fasta(fa):
    """ProteinMPNN .fa: first record is the native/backbone, rest are designs."""
    out, hdr, seq = [], None, []
    for line in open(fa):
        line = line.strip()
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(seq)))
            hdr, seq = line[1:], []
        elif line:
            seq.append(line)
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


def thread(backbone_pdb, seq, out_pdb):
    """Write binder chain A with `seq` on its backbone; target chains verbatim."""
    binder_res, lines_by_chain = [], []
    for l in open(backbone_pdb):
        if l[:4] != "ATOM":
            continue
        lines_by_chain.append(l)
        if l[21] == "A":
            n = int(l[22:26])
            if not binder_res or binder_res[-1] != n:
                binder_res.append(n)

    if len(binder_res) != len(seq):
        raise ValueError(f"{os.path.basename(backbone_pdb)}: backbone has "
                         f"{len(binder_res)} binder residues, sequence has {len(seq)}")
    aamap = {n: ONE2THREE[a] for n, a in zip(binder_res, seq)}

    serial = 0
    with open(out_pdb, "w") as f:
        for l in lines_by_chain:
            if l[21] == "A":
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
    p.add_argument("--pdbdir", required=True, help="RFdiffusion output dir")
    p.add_argument("--out", required=True)
    p.add_argument("--nseq", type=int, default=8)
    p.add_argument("--temp", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=37)
    p.add_argument("--omit", default="XC",
                   help="amino acids to omit; C excluded by default (free "
                        "cysteines in a secreted de novo binder invite "
                        "mispairing and aggregation)")
    args = p.parse_args()

    pdbs = sorted(glob.glob(os.path.join(args.pdbdir, "*.pdb")))
    if not pdbs:
        sys.exit(f"no PDBs in {args.pdbdir}")
    print(f"[mpnn] {len(pdbs)} backbones, {args.nseq} seqs each "
          f"-> {len(pdbs)*args.nseq} designs")

    t_mpnn = run_mpnn(args.pdbdir, args.out, args.nseq, args.temp,
                      args.seed, args.omit)
    print(f"[mpnn] ProteinMPNN {t_mpnn:.1f}s "
          f"({t_mpnn/max(1,len(pdbs)*args.nseq):.2f}s/design)")

    threaded = os.path.join(args.out, "threaded")
    os.makedirs(threaded, exist_ok=True)

    rows, n_ok, n_fail = [], 0, 0
    t0 = timer()
    for fa in sorted(glob.glob(os.path.join(args.out, "seqs", "*.fa"))):
        tag = os.path.basename(fa)[:-3]
        bb = os.path.join(args.pdbdir, tag + ".pdb")
        recs = parse_fasta(fa)
        for i, (hdr, seq) in enumerate(recs):
            if i == 0:
                continue  # record 0 is the input backbone, not a design
            name = f"{tag}_s{i}"
            try:
                thread(bb, seq, os.path.join(threaded, name + ".pdb"))
                n_ok += 1
            except Exception as e:
                print(f"[thread] {name} FAILED: {e}")
                n_fail += 1
                continue
            fields = dict(kv.split("=", 1) for kv in hdr.split(", ")
                          if "=" in kv)
            rows.append({
                "design": name, "parent_backbone": tag, "sequence": seq,
                "length": len(seq),
                "mpnn_score": fields.get("score", ""),
                "mpnn_seq_recovery": fields.get("seq_recovery", ""),
                "global_score": fields.get("global_score", ""),
            })
    t_thread = timer() - t0

    man = os.path.join(args.out, "mpnn_manifest.csv")
    with open(man, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["design", "parent_backbone", "sequence", "length",
                            "mpnn_score", "mpnn_seq_recovery", "global_score"])
        w.writeheader()
        w.writerows(rows)

    print(f"[thread] {n_ok} written, {n_fail} failed, {t_thread:.1f}s")
    print(f"[out] {threaded}")
    print(f"[out] {man}")


if __name__ == "__main__":
    main()
