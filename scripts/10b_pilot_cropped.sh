#!/usr/bin/env bash
# Pilot v2: same as 10_pilot_rfdiff.sh but against the 15 A cropped targets.
# comp-02 measured 1.84 -> 1.09 min/backbone from cropping; this checks the
# equivalent win here and whether interface quality holds up.
set -uo pipefail
cd "$(dirname "$0")/.."
N="${N:-8}"
OUT="pilot/rfdiff_crop"; mkdir -p "$OUT" logs
[ -f "$OUT/timings.csv" ] || echo "tag,wall_s,n_designs,s_per_design" > "$OUT/timings.csv"
run(){
  local tag="$1" pdb="$2" contigs="$3" hot="$4"
  echo "=== $tag : $N designs (cropped) ==="
  local t0=$(date +%s)
  ./scripts/rfdiffusion.sh \
    inference.output_prefix="/work/$OUT/$tag" \
    inference.input_pdb="/work/$pdb" \
    "contigmap.contigs=$contigs" "ppi.hotspot_res=$hot" \
    inference.num_designs="$N" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    hydra.run.dir=/work/.hydra > "logs/rfdiffcrop_$tag.log" 2>&1
  local rc=$? t1=$(date +%s) n=$(ls "$OUT/${tag}"_*.pdb 2>/dev/null|wc -l)
  echo "$tag rc=$rc wall=$((t1-t0))s n=$n per_design=$(( n>0 ? (t1-t0)/n : 0 ))s"
  echo "$tag,$((t1-t0)),$n,$(( n>0 ? (t1-t0)/n : 0 ))" >> "$OUT/timings.csv"
}
run classC target/glp1r_classC_crop.pdb "$(cat target/contigs_classC.txt)" '[A32,A68,A39,A69,A35]'
run classA target/glp1r_classA_crop.pdb "$(cat target/contigs_classA.txt)" '[A299,A141,A384,A298,A201]'
echo; echo "=== timings ==="; column -t -s, "$OUT/timings.csv"
