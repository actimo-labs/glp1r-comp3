# Run log 00 — Environment audit (brief §0)

Recorded 2026-10-25 ~04:12–04:20 UTC. Host: Linux 5.10.0-46-cloud-amd64.
Nothing in this log is assumed from the brief; every line was verified by command.

## Directory layout

`~/tnf-track3` is **not** a git repository (`fatal: not a git repository`).
`~/tnf-track3/glp1r-comp3` **is** its own repository → per brief §0, work on `main`.
No branch `glp1r-challenge03` needed.

```
~/tnf-track3/                      (plain dir, 11 GB, challenge-02 TNF-α artifacts)
├── RFdiffusion/                   git checkout + 3 local ckpts
├── ProteinMPNN/                   git checkout + full weight set
├── boltz/ af2/ mpnn/ rfdiff*/     comp-02 outputs (af2/ is EMPTY)
├── scripts/01..07_*.py            comp-02 pipeline, reusable code
├── target/ tnf_AB.pdb             challenge-02 target — DO NOT TOUCH
├── venv-phscan/                   biopython 1.88 + propka 3.5.1
├── rfdiffusion.sh                 Docker wrapper
└── glp1r-comp3/                   ← this repo (git, remote actimo-labs/glp1r-comp3)
```

`glp1r-comp3` before this session held only: README.md, LICENSE, CLAUDE.md,
.gitignore, glp1r-handoff-brief.md, .claude/settings.json. HEAD `d47ced0`,
branch `main`, clean, remote `https://github.com/actimo-labs/glp1r-comp3.git`.
Last 4 commits: d47ced0, 75a2f41, 550a1f0, 81cff4b.

## Tool checkouts and commits

| Tool | Path | Commit | Date | Dirty |
|---|---|---|---|---|
| RFdiffusion (checkout) | `~/tnf-track3/RFdiffusion` | `86507b6538f51fce57b5a72477165f03999ed7ae` | 2026-07-15 | clean |
| RFdiffusion (**what runs**) | docker `rosettacommons/rfdiffusion:latest` (`33d504109509`, 12.2 GB) | `fa340147b9006156b251d1ad0391e3ea8e5f73eb` | 2025-06-22 | — |
| ProteinMPNN | `~/tnf-track3/ProteinMPNN` | `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57` | 2023-06-27 | clean |
| Boltz | pip pkg in conda env `boltz` | v2.2.1 | — | — |

The host checkout of RFdiffusion is **not** the code that executes. `rfdiffusion.sh`
runs the Docker image, which is pinned a year behind the checkout. Record the
*image* commit in methods, not the checkout commit.

## Model weights on disk

| Weights | Location | Status |
|---|---|---|
| RFdiffusion Base, Complex_base, Complex_beta | `~/tnf-track3/RFdiffusion/models/` (1.4 GB) | present |
| RFdiffusion ActiveSite, Base_epoch8, Complex_Fold_base, InpaintSeq_Fold | **inside the Docker image** `/app/RFdiffusion/models` | present, but masked — see finding F1 |
| ProteinMPNN vanilla / soluble / ca (v_48_002/010/020/030) | `ProteinMPNN/*_model_weights/` | complete |
| Boltz-2 `boltz2_conf.ckpt` (2.29 GB), `boltz2_aff.ckpt` (2.06 GB), mols | `~/.boltz` (7.6 GB) | present |
| AlphaFold2 / ColabFold params | — | **ABSENT** (`af2/` empty, no `~/params`, no `~/.cache/colabfold`) |

## Environments

conda is at `~/miniforge3` (not on PATH in a fresh shell; `source env.sh` fixes it).

| Env | Python | Key packages |
|---|---|---|
| `mpnn` | 3.10.22 | torch 2.5.1+cu121 (**CUDA OK**, sees A100), biopython 1.88, numpy 2.2.6 |
| `boltz` | 3.11.17 | boltz 2.2.1, torch 2.14.1+cu126 (**CUDA OK**), rdkit 2026.09.1, biopython 1.84 |
| `viz` | — | not inspected |
| `venv-phscan` | 3.10 | biopython 1.88, propka 3.5.1, scipy 1.15.3 |

Absent from every env: **pyrosetta, MDAnalysis, ProDy, jax/alphafold, openmm, pdbfixer**.
Absent binaries: `colabfold_batch`, `foldseek`, `mmseqs`, `pymol`. `tmux` and `docker` present
(user is in the `docker` group; GPU passthrough verified inside the container).

## GPU

`NVIDIA A100-SXM4-80GB`, driver 550.90.07, CUDA 12.4. **Completely idle** —
1 MiB / 81920 MiB used, 0% util, no processes, 34 °C. Exclusive for this run.
Verified usable from inside the RFdiffusion container (`--gpus all` → GPU 0 listed).

## Disk — the binding constraint

```
/dev/sda1  2.0T  1.9T used  37G avail  99%
```

The volume is shared. **1.8 TB belongs to another user** (`/home/jaberjaber`:
boltzgen training targets 42 GB, MD trajectory 23 GB, HF caches, etc.). Not ours
to delete. `~/HarjasG` totals ~29 GB (18 GB miniforge3 + 11 GB tnf-track3).

Working headroom for this campaign: **37 GB, and it can shrink without warning**
if the other user writes. Challenge-02 consumed ~8.3 GB of outputs
(rfdiff 4.3 GB + boltz 4.0 GB) for a comparable campaign.

## Network

All reachable, HTTP 200:

| Endpoint | Time |
|---|---|
| `files.rcsb.org/download/6X18.pdb` | 0.57 s |
| `www.rcsb.org` | 0.83 s |
| `proteinbase.com/.../challenges/glp1r` | 1.41 s |
| `github.com` | 0.04 s |
| `pypi.org/simple/` | 0.11 s |

## 6X18 verified against the brief

Downloaded (974 KB). **Matches the brief.** ELECTRON MICROSCOPY, 2.10 Å.
Chains: `A` Gαs-short, `B` Gβ1, `G` Gγ2, `N` Nanobody35, `P` GLP-1(7-36)NH2,
`R` GLP-1 receptor. So chain R is the receptor as stated, chain P is the native
ligand available for deriving hotspots by contact (brief §6: compute, don't hardcode),
and Gs+Nb35 are present if we choose to retain them to hold the active intracellular face.

## Git housekeeping done

- `user.name` / `user.email` were **unset globally** → commits would have failed.
  Set repo-locally (`HarjasG` / `HarjasG@users.noreply.github.com`). Trivial fix, per brief §0.
- `.gitignore` is the stock Python template and ignored `target/` (Python packaging
  artifact dir) and `*.log`. Both collide with this project. Appended negations plus
  rules to keep large artifacts out of git.
- Created `runlog/ target/ scripts/ manifests/`.
