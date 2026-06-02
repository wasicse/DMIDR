#!/usr/bin/env bash
# Apptainer wrapper: runs GROMACS MD simulation inside gmx-gpu.sif.
# Usage: bash run_gromacs.sh <wt.pdb> <mut.pdb> <outdir>
#
# Required env:
#   GMX_SIF   — absolute path to gmx-gpu.sif
#
# Optional env (passed into container):
#   GMX_STEPS        (default 3000000)
#   GMX_SAVE_EVERY   (default 1000)
#   GMX_FF           (default amber14sb)
#   GMX_WATER        (default tip3p)
#   GMX_TEMP         (default 280)
#   GMX_SALT         (default 0.15)
#   GMX_PADDING      (default 2.5)
#   GMX_GPU_OPTS     (default "-nb gpu -pme gpu -bonded gpu -pin on -ntmpi 1 -ntomp 10")
set -euo pipefail

ORIG="$(realpath "${1:?arg1: path to WT PDB}")"
MUT="$(realpath  "${2:?arg2: path to mutant PDB}")"
OUT="$(realpath  "${3:?arg3: output directory}")"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MDP_DIR="$PROJECT_ROOT/configs/gromacs"

GMX_SIF="${GMX_SIF:?GMX_SIF must be set — path to gmx-gpu.sif}"

mkdir -p "$OUT"

apptainer exec --nv \
  -B "$PROJECT_ROOT:$PROJECT_ROOT" \
  --env STEPS="${GMX_STEPS:-3000000}" \
  --env SAVE_EVERY="${GMX_SAVE_EVERY:-1000}" \
  --env FF="${GMX_FF:-amber14sb}" \
  --env WATER="${GMX_WATER:-tip3p}" \
  --env TEMP="${GMX_TEMP:-280}" \
  --env SALT="${GMX_SALT:-0.15}" \
  --env PADDING="${GMX_PADDING:-2.5}" \
  --env MD_GPU_OPTS="${GMX_GPU_OPTS:--nb gpu -pme gpu -pin on -ntmpi 1 -ntomp $(nproc)}" \
  "$GMX_SIF" \
  bash "$SCRIPT_DIR/run_gromacs_inner.sh" "$ORIG" "$MUT" "$OUT" "$MDP_DIR"
