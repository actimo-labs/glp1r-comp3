#!/usr/bin/env bash
# Triage: generated backbones -> ranked candidates ready for Boltz-2.
#
#   interface QC  ->  ProteinMPNN  ->  thread onto FULL receptor
#                 ->  de novo check  ->  AF2-ig (class C only)
#
# Per-class gating, from runlog/04-af2-triage-calibration.md:
#
#   class C  targets the soluble ECD, where AF2 is reliable. AF2-ig IS a gate:
#            pae_interaction < 22.0 (below all 8 negatives; the native GLP-1
#            positive scores 17.39 in this backbone-only representation).
#
#   class A  target the TMD, where AF2-ig has median DockQ 0.03 on GPCR
#   class B  complexes (brief section 3). AF2-ig is NOT a gate for these -- its
#            pae_interaction saturates at ~27.5 for every design, which is "no
#            information", not "no interface". They are pre-filtered on
#            interface QC + MPNN score + plddt_binder (does the binder fold at
#            all) and scored by Boltz-2 instead.
set -uo pipefail
cd "$(dirname "$0")/.."

GEN="${GEN:-gen/rfdiff}"
OUT="${OUT:-gen/triage}"
NSEQ="${NSEQ:-4}"
mkdir -p "$OUT" logs

echo "### 1. interface QC"
# hotspots per class, original receptor numbering
qc(){ # tag hotspots target_start
  local tag="$1"
  /home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/12_interface_qc.py \
    --dir "$GEN" --tag "$tag" --hotspots "$2" --target-start "$3" \
    2>&1 | tail -6
}
qc classC 32,68,39,69,35       29
qc classA 299,141,384,298,201 136
qc classB 32,68,39,69,35,299,141 29

echo
echo "### 2. collect QC survivors"
/home/HarjasG/miniforge3/envs/af2ig/bin/python - "$GEN" "$OUT" <<'PY'
import csv, glob, os, shutil, sys
gen, out = sys.argv[1], sys.argv[2]
keep_dir = os.path.join(out, "qc_pass"); os.makedirs(keep_dir, exist_ok=True)
n_in = n_keep = 0
for f in glob.glob(os.path.join(gen, "interface_qc_*.csv")):
    for r in csv.DictReader(open(f)):
        n_in += 1
        if not r["verdict"].startswith("ok"):
            continue
        src = os.path.join(gen, r["design"] + ".pdb")
        if os.path.isfile(src):
            shutil.copy(src, keep_dir)
            trb = src[:-4] + ".trb"
            if os.path.isfile(trb):
                shutil.copy(trb, keep_dir)
            n_keep += 1
print(f"  {n_keep}/{n_in} backbones pass QC -> {keep_dir}")
PY

echo
echo "### 3. ProteinMPNN + thread onto the FULL receptor"
/home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/11_mpnn_thread.py \
  --pdbdir "$OUT/qc_pass" --out "$OUT/mpnn" --nseq "$NSEQ" --full-target \
  2>&1 | tail -6

echo
echo "### 4. de novo compliance"
/home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/20_denovo_check.py \
  --manifest "$OUT/mpnn/mpnn_manifest.csv" \
  --out "$OUT/denovo_check.csv" 2>&1 | tail -6 || true

echo
echo "### 5. AF2-ig on class C only"
mkdir -p "$OUT/af2_in"
cp "$OUT"/mpnn/threaded/classC*.pdb "$OUT/af2_in/" 2>/dev/null
n=$(ls "$OUT/af2_in"/*.pdb 2>/dev/null | wc -l)
echo "  $n class C designs"
if [ "$n" -gt 0 ]; then
  /home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/af2ig/predict_pdb.py \
    -pdbdir "$OUT/af2_in" \
    -outpdbdir "$OUT/af2_out" \
    -scorefilename "$OUT/af2.sc" \
    -checkpoint_name "$OUT/af2.check" > logs/triage_af2.log 2>&1
  echo "  af2 rc=$?  $(grep -c '^SCORE:' "$OUT/af2.sc") lines"
fi

echo
echo "### done -> $OUT"
echo "next: scripts/41_select.py to build the Boltz shortlist"
