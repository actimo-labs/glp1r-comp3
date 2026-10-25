#!/usr/bin/env bash
# Pilot RFdiffusion batches for class A (vestibule rim) and class C (ECD groove).
# Records wall-clock per design per class. Run from the repo root.
#
# Class B (two-domain mimetic) is NOT here -- it needs motif scaffolding of the
# message module and is set up separately once A and C have validated the stack.
set -uo pipefail
cd "$(dirname "$0")/.."

N="${N:-8}"
OUT="pilot/rfdiff"
mkdir -p "$OUT" logs

run() {
  local tag="$1" pdb="$2" contigs="$3" hotspots="$4"
  echo "=== $tag : $N designs ==="
  local t0=$(date +%s)
  ./scripts/rfdiffusion.sh \
    inference.output_prefix="/work/$OUT/$tag" \
    inference.input_pdb="/work/$pdb" \
    "contigmap.contigs=$contigs" \
    "ppi.hotspot_res=$hotspots" \
    inference.num_designs="$N" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    hydra.run.dir=/work/.hydra > "logs/rfdiff_$tag.log" 2>&1
  local rc=$? t1=$(date +%s)
  local n=$(ls "$OUT/${tag}"_*.pdb 2>/dev/null | wc -l)
  echo "$tag rc=$rc  wall=$((t1-t0))s  designs=$n  per_design=$(( n>0 ? (t1-t0)/n : 0 ))s"
  echo "$tag,$((t1-t0)),$n,$(( n>0 ? (t1-t0)/n : 0 ))" >> "$OUT/timings.csv"
}

[ -f "$OUT/timings.csv" ] || echo "tag,wall_s,n_designs,s_per_design" > "$OUT/timings.csv"

# ---- class C : ECD groove. Target 29-129, contiguous. --------------------
run classC target/glp1r_ECD.pdb \
    '[A29-129/0 60-90]' \
    '[A32,A68,A39,A69,A35]'

# ---- class A : vestibule rim. Target 136-423 with a gap at 339-343. ------
run classA target/glp1r_TMD.pdb \
    '[A136-338/A344-423/0 60-90]' \
    '[A299,A141,A384,A298,A201]'

echo
echo "=== timings ==="
column -t -s, "$OUT/timings.csv"
