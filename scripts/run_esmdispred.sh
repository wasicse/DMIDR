#!/usr/bin/env bash
# Run ESMDisPred (Docker) for each mutant FASTA file in FASTA_DIR.
# One .caid file is produced per mutation block size.
#
# Usage:
#   bash scripts/run_esmdispred.sh \
#       <fasta_dir> <output_dir> <image> <large_models_dir> <seq_name> \
#       [max_mutations=5] [model=3]
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: bash scripts/run_esmdispred.sh <fasta_dir> <output_dir> <image> <large_models_dir> <seq_name> [max_mutations] [model]" >&2
  exit 1
fi

FASTA_DIR="$1"
OUTPUT_DIR="$2"
IMAGE_NAME="$3"
LARGE_MODELS_DIR="$(realpath "$4")"
SEQ_NAME="$5"
MAX_MUTATIONS="${6:-5}"
MODEL="${7:-3}"

ESMDISPRED_DIR="$(realpath "$(dirname "$0")/tools/ESMDisPred")"
ESMpath="/opt/ESMDisPred"

if [[ ! -d "$ESMDISPRED_DIR" ]]; then
  echo "Error: ESMDisPred tool not found at $ESMDISPRED_DIR" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

for i in $(seq 1 "$MAX_MUTATIONS"); do
  file=$(find "$FASTA_DIR" -maxdepth 1 -type f -name "*_mutants_${i}res.fasta" | head -n 1)
  if [[ -z "$file" ]]; then
    echo "Skipping mutation size ${i}: no FASTA file found"
    continue
  fi
  file="$(realpath "$file")"
  fasta_filename="$(basename "$file")"

  output_name="${SEQ_NAME}_mutants_${i}res.caid"
  tmp_out="$(mktemp -d)"

  pushd "$ESMDISPRED_DIR" > /dev/null
  mkdir -p features

  docker run --rm \
    --gpus all \
    --user "$(id -u):$(id -g)" \
    -e HOME="$ESMpath" \
    -e XDG_CACHE_HOME="$ESMpath/.cache" \
    -e TORCH_HOME="$ESMpath/largeModels" \
    -v "$file":"$ESMpath/example/$fasta_filename":ro \
    -v "$tmp_out":"$ESMpath/outputs":rw \
    -v "$(pwd)/features":"$ESMpath/features":rw \
    -v "$LARGE_MODELS_DIR":"$ESMpath/largeModels":rw \
    -v "$(pwd)/run_ESMDisPred.sh":"$ESMpath/run_ESMDisPred.sh":ro \
    -v "$(pwd)/scripts/run_Dispredict3.sh":"$ESMpath/scripts/run_Dispredict3.sh":ro \
    -v "$(pwd)/scripts/run_ESMDisPred.py":"$ESMpath/scripts/run_ESMDisPred.py":ro \
    -v "$(pwd)/scripts/run_ESM2.py":"$ESMpath/scripts/run_ESM2.py":ro \
    -v "$(pwd)/tools/Dispredict3.0/tools/fldpnn/run_flDPnn.py":"$ESMpath/tools/Dispredict3.0/tools/fldpnn/run_flDPnn.py":ro \
    -v "$(pwd)/scripts/transformer_Inference.py":"$ESMpath/scripts/transformer_Inference.py":ro \
    -v "$(pwd)/scripts/preprocess.py":"$ESMpath/scripts/preprocess.py":ro \
    -v "$(pwd)/models":"$ESMpath/models":ro \
    -v "$(pwd)/requirements.txt":"$ESMpath/requirements.txt":ro \
    -v "$(pwd)/run_downloadLargeModels.sh":"$ESMpath/run_downloadLargeModels.sh":ro \
    "$IMAGE_NAME" \
    ./run_ESMDisPred.sh "$ESMpath/example/$fasta_filename" outputs "$MODEL"

  popd > /dev/null

  mapfile -t caid_files < <(find "$tmp_out" -name "*.caid" | sort)
  if [[ ${#caid_files[@]} -eq 0 ]]; then
    echo "ERROR: No .caid output for mutation size ${i}" >&2
    rm -rf "$tmp_out"
    exit 1
  fi
  cat "${caid_files[@]}" > "${OUTPUT_DIR}/${output_name}"
  rm -rf "$tmp_out"
  echo "Created ${OUTPUT_DIR}/${output_name}"
done
