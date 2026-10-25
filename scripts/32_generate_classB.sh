#!/usr/bin/env bash
# Class B generation -- COMPLIANT formulation.
#
# WHY THIS CHANGED
# ----------------
# The first version held GLP-1 residues 7-15 BACKBONE fixed as an RFdiffusion
# motif, with the sequence designed de novo. Competition FAQ 8:
#
#   "They must be de novo and zero-shot, defined as not having relied on an
#    initial binder as a starting point ... having adequate sequence- AND
#    structural-diversity [from known proteins] ... You may not take an existing
#    binder and modify it. Designs must be produced from scratch."
#
# Known binders are explicitly permitted "to evaluate/calibrate design filtering
# metrics" -- which is exactly how the negative controls use comp-02's TNF
# designs, and that stays. They are NOT permitted as a starting point.
#
# GLP-1 is the native binder of this exact target. A 9-residue backbone graft at
# the functional site has full sequence diversity but deliberately REDUCED
# structural diversity, precisely where it matters. That fails the
# structural-diversity clause on a strict reading, and 10 of 20 slots is too
# much to stake on a permissive one.
#
# COMPLIANT APPROACH
# ------------------
# No GLP-1 coordinates anywhere. The receptor's own pocket defines the target:
# RFdiffusion is pointed at the deep activation-trigger residues via
# ppi.hotspot_res and must discover its own way in. Binders are longer than
# class A/C so there is enough chain to span the vestibule.
#
#   A387 (35.1 A)  the free-alpha-amine partner, the activation trigger
#   A234 (29.3 A)  direct peptide contact
#   A233 (29.0 A)  deep pocket
#   A391 (38.4 A)  deep
#   A241 (39.5 A)  deepest
#
# Everything here is derived from the TARGET structure (6X18 chain R), which the
# competition lists as the structure reference. Nothing is taken from the ligand.
set -uo pipefail
cd "$(dirname "$0")/.."

NB="${NB:-250}"
OUT="${OUT:-gen/rfdiff}"
LENS=(${LENS_B:-90 110 130})
mkdir -p "$OUT" logs
[ -f "$OUT/timings.csv" ] || echo "tag,len,wall_s,n_designs,s_per_design" > "$OUT/timings.csv"

BSEG="$(sed -E 's/^\[//; s!/0 .*$!!' target/contigs_classB2_cropONLY.txt)"
HOT='[A387,A234,A233,A391,A241]'
echo "class B compliant: target=$BSEG"
echo "                   hotspots=$HOT (deep pocket, receptor-derived only)"

for L in "${LENS[@]}"; do
  n=$(( NB / ${#LENS[@]} ))
  [ "$n" -le 0 ] && continue
  echo "=== classB len=$L : $n designs ==="
  t0=$(date +%s)
  ./scripts/rfdiffusion.sh \
    inference.output_prefix="/work/$OUT/classBL${L}" \
    inference.input_pdb=/work/target/glp1r_classB2_crop.pdb \
    "contigmap.contigs=[${BSEG}/0 ${L}-${L}]" \
    "ppi.hotspot_res=$HOT" \
    inference.num_designs="$n" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    inference.write_trajectory=False \
    hydra.run.dir=/work/.hydra >> "logs/gen_classBL${L}.log" 2>&1
  rc=$?; t1=$(date +%s)
  got=$(ls "$OUT/classBL${L}"_*.pdb 2>/dev/null | wc -l)
  echo "classB len=$L rc=$rc wall=$((t1-t0))s designs=$got per_design=$(( got>0 ? (t1-t0)/got : 0 ))s"
  echo "classB,$L,$((t1-t0)),$got,$(( got>0 ? (t1-t0)/got : 0 ))" >> "$OUT/timings.csv"
done
echo "CLASSB_DONE"
