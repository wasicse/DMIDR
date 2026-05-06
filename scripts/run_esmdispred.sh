#!/usr/bin/env bash
# Run ESMDisPred (Docker) for each mutant FASTA file in FASTA_DIR.
# One .caid file is produced per mutation block size, concatenating all
# per-sequence .caid outputs from ESMDisPred.
#
# Usage:
#   bash scripts/run_esmdispred.sh \
#       <fasta_dir> <output_dir> <image> <large_models_dir> <seq_name> \
#       [max_mutations=5] [model=3]
#
# model argument: 1=ESMDisPred-1  2=ESMDisPred-2  3=ESMDisPred-2PDB
#                 4=ESMDisPred-DNN  all=run all variants
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: bash scripts/run_esmdispred.sh <fasta_dir> <output_dir> <image> <large_models_dir> <seq_name> [max_mutations] [model]" >&2
  exit 1
fi

FASTA_DIR="$1"
OUTPUT_DIR="$2"
IMAGE_NAME="$3"
LARGE_MODELS_DIR="$4"
SEQ_NAME="$5"
MAX_MUTATIONS="${6:-5}"
MODEL="${7:-3}"

mkdir -p "$OUTPUT_DIR"

for i in $(seq 1 "$MAX_MUTATIONS"); do
  file=$(find "$FASTA_DIR" -maxdepth 1 -type f -name "*_mutants_${i}res.fasta" | head -n 1)
  if [[ -z "$file" ]]; then
    echo "Skipping mutation size ${i}: no FASTA file found"
    continue
  fi

  output_name="${SEQ_NAME}_mutants_${i}res.caid"
  tmp_out="${OUTPUT_DIR}/.esmdispred_tmp_${i}"
  mkdir -p "$tmp_out"

  docker run --rm \
    -v "$(realpath "$file"):/opt/ESMDisPred/example/sample.fasta:ro" \
    -v "$(realpath "$LARGE_MODELS_DIR"):/opt/ESMDisPred/largeModels:ro" \
    -v "$(realpath "$tmp_out"):/opt/ESMDisPred/outputs:rw" \
    "$IMAGE_NAME" \
    bash -lc "/opt/ESMDisPred/run_ESMDisPred.sh \
      /opt/ESMDisPred/example/sample.fasta \
      /opt/ESMDisPred/outputs \
      ${MODEL}"

  # Concatenate all per-sequence .caid files into one file for this block size
  mapfile -t caid_files < <(find "$tmp_out" -name "*.caid" | sort)
  if [[ ${#caid_files[@]} -eq 0 ]]; then
    echo "ERROR: No .caid output found in $tmp_out for mutation size ${i}" >&2
    rm -rf "$tmp_out"
    exit 1
  fi
  cat "${caid_files[@]}" > "${OUTPUT_DIR}/${output_name}"

  rm -rf "$tmp_out"
  echo "Created ${OUTPUT_DIR}/${output_name}"
done
