#!/usr/bin/env bash
# Main RFdiffusion generation. Counts are parameters because the GPU budget is
# the binding constraint -- see runlog/03-pilot-and-budget.md for the arithmetic.
#
#   NA / NB / NC   backbones per class
#
# Checkpointing: RFdiffusion skips designs whose output PDB already exists, so
# re-running this script resumes rather than restarts. With ~30 h and a single
# A100 there is no budget for a lost run (brief section 0).
#
# Portfolio after the depth analysis (runlog/02-epitope-depth.md): 4A / 10B / 6C.
# Class B gets the most slots and the most backbones because it is the only
# architecture whose extended N-terminus reaches E387 at 35 A.
set -uo pipefail
cd "$(dirname "$0")/.."

NA="${NA:-150}"
NB="${NB:-250}"
NC="${NC:-200}"
OUT="${OUT:-gen/rfdiff}"
mkdir -p "$OUT" logs
[ -f "$OUT/timings.csv" ] || echo "tag,wall_s,n_designs,s_per_design" > "$OUT/timings.csv"

run(){
  local tag="$1" pdb="$2" contigs="$3" hot="$4" n="$5"
  [ "$n" -eq 0 ] && { echo "skip $tag (n=0)"; return; }
  echo "=== $tag : $n designs ==="
  local t0=$(date +%s)
  ./scripts/rfdiffusion.sh \
    inference.output_prefix="/work/$OUT/$tag" \
    inference.input_pdb="/work/$pdb" \
    "contigmap.contigs=$contigs" \
    "ppi.hotspot_res=$hot" \
    inference.num_designs="$n" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    hydra.run.dir=/work/.hydra >> "logs/gen_$tag.log" 2>&1
  local rc=$? t1=$(date +%s)
  local got=$(ls "$OUT/${tag}"_*.pdb 2>/dev/null | wc -l)
  echo "$tag rc=$rc wall=$((t1-t0))s total_designs=$got per_design=$(( got>0 ? (t1-t0)/got : 0 ))s"
  echo "$tag,$((t1-t0)),$got,$(( got>0 ? (t1-t0)/got : 0 ))" >> "$OUT/timings.csv"
}

# ---- class C : ECD groove. Cheapest per design, published precedent. -------
run classC target/glp1r_classC_crop.pdb \
    "$(cat target/contigs_classC.txt)" '[A32,A68,A39,A69,A35]' "$NC"

# ---- class A : vestibule rim. Thin 6-residue epitope, insurance only. ------
run classA target/glp1r_classA_crop.pdb \
    "$(cat target/contigs_classA.txt)" '[A299,A141,A384,A298,A201]' "$NA"

# ---- class B : two-domain mimetic. Motif-scaffolded. -----------------------
# The GLP-1 7-15 BACKBONE is held fixed in its native frame as chain P, so the
# designed N-terminus starts at native insertion depth with a free alpha-amine
# pointing at E387. Sequence is designed from scratch downstream; no GLP-1
# sequence survives, which scripts/20_denovo_check.py verifies.
# Contig: target segments, chain break, then motif + 40-80 designed residues.
run classB target/glp1r_classB2_crop_motif.pdb \
    '[A29-53/A62-129/A136-155/A187-241/A280-317/A361-397/0 P7-15/40-80]' \
    '[A32,A68,A39,A69,A35,A299,A141]' "$NB"

echo
echo "=== timings ==="
column -t -s, "$OUT/timings.csv"
echo "total backbones: $(ls "$OUT"/*.pdb 2>/dev/null | wc -l)"
