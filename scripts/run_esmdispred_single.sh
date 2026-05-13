#!/usr/bin/env bash
# Run ESMDisPred on a single FASTA file and write one .caid output file.
#
# Usage:
#   bash scripts/run_esmdispred_single.sh \
#     <fasta_file> <output_dir> <image> <large_models_dir> <output_stem> [model=3]
#
# model: 1=ESMDisPred-1  2=ESMDisPred-2  3=ESMDisPred-2PDB (default)
#        4=ESMDisPred-DNN  all=run all variants
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <fasta_file> <output_dir> <image> <large_models_dir> <output_stem> [model]" >&2
  exit 1
fi

FASTA="$1"
OUTPUT_DIR="$2"
IMAGE_NAME="$3"
LARGE_MODELS_DIR="$4"
OUTPUT_STEM="$5"
MODEL="${6:-3}"

mkdir -p "$OUTPUT_DIR"
tmp_out="$OUTPUT_DIR/.esmdispred_single_$$"
mkdir -p "$tmp_out"

docker run --rm \
  -v "$(realpath "$FASTA"):/opt/ESMDisPred/example/sample.fasta:ro" \
  -v "$(realpath "$LARGE_MODELS_DIR"):/opt/ESMDisPred/largeModels:ro" \
  -v "$(realpath "$tmp_out"):/opt/ESMDisPred/outputs:rw" \
  "$IMAGE_NAME" \
  bash -lc "/opt/ESMDisPred/run_ESMDisPred.sh \
    /opt/ESMDisPred/example/sample.fasta \
    /opt/ESMDisPred/outputs \
    ${MODEL}"

mapfile -t caid_files < <(find "$tmp_out" -name "*.caid" | sort)
if [[ ${#caid_files[@]} -eq 0 ]]; then
  echo "ERROR: No .caid output from ESMDisPred for $FASTA" >&2
  rm -rf "$tmp_out"
  exit 1
fi

cat "${caid_files[@]}" > "${OUTPUT_DIR}/${OUTPUT_STEM}.caid"
rm -rf "$tmp_out"
echo "Created ${OUTPUT_DIR}/${OUTPUT_STEM}.caid"
