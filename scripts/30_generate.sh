#!/usr/bin/env bash
# Main RFdiffusion generation.
#
#   NA / NB / NC   backbones per class (split evenly across that class's lengths)
#
# WHY DISCRETE LENGTHS INSTEAD OF A RANGE
# ---------------------------------------
# The obvious contig is `60-90`, letting RFdiffusion pick a length per design.
# That is expensive downstream: JAX retraces AF2 for every new input shape.
# Measured on this box: 122 s on a new complex length, 10.5 s on a repeat. A
# continuous 60-90 range gives ~31 distinct lengths per class, so ~31 compiles
# per class -- about an hour of pure overhead each, before any inference.
#
# Sampling 3 discrete lengths per class instead gives 3 compiles per class
# (~6 min each) and costs nothing in design quality: length diversity is
# quantised, not removed, and the brief only asks for >=3 backbone topologies.
#
# Checkpointing: RFdiffusion skips designs whose output PDB already exists, so
# re-running resumes rather than restarts.
#
# Portfolio after the depth analysis (runlog/02-epitope-depth.md): 4A / 10B / 6C.
set -uo pipefail
cd "$(dirname "$0")/.."

NA="${NA:-120}"
NB="${NB:-250}"
NC="${NC:-200}"
OUT="${OUT:-gen/rfdiff}"
mkdir -p "$OUT" logs
[ -f "$OUT/timings.csv" ] || echo "tag,len,wall_s,n_designs,s_per_design" > "$OUT/timings.csv"

# class A and C: plain binder lengths
LENS_AC=(65 75 85)
# class B: 9-residue motif is prepended, so these are the DESIGNED residues.
# Binder totals come out 54 / 69 / 84.
LENS_B=(45 60 75)

run(){
  local tag="$1" pdb="$2" contigs="$3" hot="$4" n="$5" len="$6"
  [ "$n" -le 0 ] && return
  echo "=== $tag len=$len : $n designs ==="
  local t0=$(date +%s)
  ./scripts/rfdiffusion.sh \
    inference.output_prefix="/work/$OUT/${tag}L${len}" \
    inference.input_pdb="/work/$pdb" \
    "contigmap.contigs=$contigs" \
    "ppi.hotspot_res=$hot" \
    inference.num_designs="$n" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    hydra.run.dir=/work/.hydra >> "logs/gen_${tag}L${len}.log" 2>&1
  local rc=$? t1=$(date +%s)
  local got=$(ls "$OUT/${tag}L${len}"_*.pdb 2>/dev/null | wc -l)
  echo "$tag len=$len rc=$rc wall=$((t1-t0))s designs=$got per_design=$(( got>0 ? (t1-t0)/got : 0 ))s"
  echo "$tag,$len,$((t1-t0)),$got,$(( got>0 ? (t1-t0)/got : 0 ))" >> "$OUT/timings.csv"
}

# Segments are DERIVED from the crop outputs, never retyped. Hardcoding them
# once already produced an off-by-one (A62-74 for a crop that says A63-74),
# which RFdiffusion would have accepted silently as a different epitope.
seg_of(){ sed -E 's/^\[//; s!/0 .*$!!' "$1"; }
CSEG="$(seg_of target/contigs_classC.txt)"
ASEG="$(seg_of target/contigs_classA.txt)"
BSEG="$(seg_of target/contigs_classB2_cropONLY.txt)"
echo "derived segments:"
echo "  C: $CSEG"
echo "  A: $ASEG"
echo "  B: $BSEG"

# ---- class C : ECD groove. Cheapest, published precedent, 6 slots. --------
for L in "${LENS_AC[@]}"; do
  run classC target/glp1r_classC_crop.pdb "[${CSEG}/0 ${L}-${L}]" \
      '[A32,A68,A39,A69,A35]' $(( NC / ${#LENS_AC[@]} )) "$L"
done

# ---- class A : vestibule rim. Thin epitope, 4 slots, insurance. -----------
for L in "${LENS_AC[@]}"; do
  run classA target/glp1r_classA_crop.pdb "[${ASEG}/0 ${L}-${L}]" \
      '[A299,A141,A384,A298,A201]' $(( NA / ${#LENS_AC[@]} )) "$L"
done

# ---- class B : two-domain mimetic, motif-scaffolded, 10 slots. ------------
# GLP-1 7-15 BACKBONE held in its native frame as chain P, so the designed
# N-terminus starts at native insertion depth with a free alpha-amine aimed at
# E387. Verified in the smoke test: motif lands at binder residues 1-9 and
# reaches 30.4 A depth (native His7 = 31.2 A). Sequence is designed from
# scratch downstream; 20_denovo_check.py verifies no GLP-1 sequence survives.
for L in "${LENS_B[@]}"; do
  run classB target/glp1r_classB2_crop_motif.pdb "[${BSEG}/0 P7-15/${L}-${L}]" \
      '[A32,A68,A39,A69,A35,A299,A141]' $(( NB / ${#LENS_B[@]} )) "$L"
done

echo
echo "=== timings ==="; column -t -s, "$OUT/timings.csv"
echo "total backbones: $(ls "$OUT"/*.pdb 2>/dev/null | wc -l)"
echo "distinct binder lengths: A/C ${LENS_AC[*]}  B $(printf '%s ' "${LENS_B[@]}") (+9 motif)"
