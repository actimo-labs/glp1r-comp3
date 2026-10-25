#!/usr/bin/env bash
# Runs the remaining pipeline unattended, so nothing waits on an interactive
# session. Each stage waits for its real prerequisite with EXACT tmux matching
# (-t =NAME); prefix matching previously deadlocked the chain.
#
#   wait for class B generation + class C AF2
#     -> interface QC on class B
#     -> ProteinMPNN on class B survivors
#     -> de novo check
#     -> select Boltz shortlist (all classes)
#     -> Boltz-2
#     -> rank + build submission package (gated, nothing uploaded)
set -uo pipefail
cd /home/HarjasG/tnf-track3/glp1r-comp3
log(){ echo "[$(date -u +%H:%M:%S)] $*"; }

log "waiting for class B generation and class C AF2"
while tmux has-session -t =genall 2>/dev/null || tmux has-session -t =af2Cecd 2>/dev/null; do
  sleep 60
done
log "prerequisites finished"

log "=== class B interface QC ==="
/home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/12_interface_qc.py \
  --dir gen/rfdiff --tag classB --hotspots 387,234,233,391,241 --target-start 29 2>&1 | tail -6

log "=== stage class B survivors ==="
/home/HarjasG/miniforge3/envs/af2ig/bin/python - <<'PY'
import csv, os, shutil
os.makedirs("gen/triage/qc_pass_B", exist_ok=True)
k=t=0
f="gen/rfdiff/interface_qc_classB.csv"
if os.path.isfile(f):
    for r in csv.DictReader(open(f)):
        t+=1
        if not r["verdict"].startswith("ok"): continue
        src=f"gen/rfdiff/{r['design']}.pdb"
        if os.path.isfile(src):
            shutil.copy(src,"gen/triage/qc_pass_B/")
            trb=src[:-4]+".trb"
            if os.path.isfile(trb): shutil.copy(trb,"gen/triage/qc_pass_B/")
            k+=1
print(f"  class B: {k}/{t} staged")
PY

log "=== ProteinMPNN on class B ==="
/home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/11_mpnn_thread.py \
  --pdbdir gen/triage/qc_pass_B --out gen/triage/mpnn_B --nseq 4 --full-target 2>&1 | tail -4

log "=== merge manifests ==="
/home/HarjasG/miniforge3/envs/af2ig/bin/python - <<'PY'
import csv, os
rows=[]
for p in ("gen/triage/mpnn/mpnn_manifest.csv","gen/triage/mpnn_B/mpnn_manifest.csv"):
    if os.path.isfile(p): rows += list(csv.DictReader(open(p)))
if rows:
    with open("gen/triage/mpnn_all.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"  {len(rows)} sequences total")
PY

log "=== de novo check (all classes) ==="
/home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/20_denovo_check.py \
  --manifest gen/triage/mpnn_all.csv --out gen/triage/denovo_all.csv 2>&1 | tail -3 || true

log "=== select Boltz shortlist ==="
cp -f gen/triage/af2_ecd.sc gen/triage/af2.sc 2>/dev/null
cp -f gen/triage/denovo_all.csv gen/triage/denovo_check.csv 2>/dev/null
/home/HarjasG/miniforge3/envs/af2ig/bin/python - <<'PY'
import csv, os, shutil
# point the selector at the merged manifest
src="gen/triage/mpnn_all.csv"; dst="gen/triage/mpnn/mpnn_manifest.csv"
if os.path.isfile(src): shutil.copy(src,dst)
PY
/home/HarjasG/miniforge3/envs/af2ig/bin/python scripts/41_select.py \
  --triage gen/triage --gen gen/rfdiff --out gen/triage \
  --n-c 60 --n-a 40 --n-b 100 2>&1 | tail -8

log "=== Boltz-2 on the shortlist ==="
/home/HarjasG/miniforge3/envs/boltz/bin/python scripts/13_boltz.py \
  --manifest gen/triage/boltz_shortlist.csv --out gen/boltz --seed 1 2>&1 | tail -12

log "=== AUTOPILOT COMPLETE ==="
log "next: rank + build submission (needs review, not automated)"
echo AUTOPILOT_DONE
