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

PROJECT_ROOT="$(realpath "$(dirname "$0")/..")"
ESMDISPRED_DIR="$(realpath "$(dirname "$0")/tools/ESMDisPred")"
DISOCOMB_SH="$PROJECT_ROOT/scripts/esmdispred/DisoComb.sh"
PSSM_CACHE_DIR="${PSSM_CACHE_DIR:-}"
ESMpath="/opt/ESMDisPred"

if [[ ! -d "$ESMDISPRED_DIR" ]]; then
  echo "Error: ESMDisPred tool not found at $ESMDISPRED_DIR" >&2
  exit 1
fi

# Build optional extra mounts for parallel PSI-BLAST and PSSM cache
extra_mounts=()
if [[ -f "$DISOCOMB_SH" ]]; then
  extra_mounts+=(
    -v "$DISOCOMB_SH":"$ESMpath/tools/Dispredict3.0/tools/fldpnn/DisoComb.sh":ro
  )
fi
if [[ -n "$PSSM_CACHE_DIR" ]]; then
  mkdir -p "$PSSM_CACHE_DIR"
  extra_mounts+=(
    -v "$PSSM_CACHE_DIR":"$ESMpath/pssm_cache":rw
    -e PSSM_CACHE="$ESMpath/pssm_cache"
  )
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
  if [[ -f "$OUTPUT_DIR/$output_name" ]]; then
    echo "  [skip] ${i}res — $output_name already exists"
    continue
  fi
  tmp_out="$(mktemp -d)"

  pushd "$ESMDISPRED_DIR" > /dev/null
  mkdir -p features

  docker run --rm \
    --gpus all \
    --user "$(id -u):$(id -g)" \
    -e HOME="$ESMpath" \
    -e XDG_CACHE_HOME="$ESMpath/.cache" \
    -e TORCH_HOME="$ESMpath/largeModels" \
    "${extra_mounts[@]}" \
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

  # Map numeric model to the subdirectory ESMDisPred writes into.
  case "$MODEL" in
    1|ESMDisPred-1)   _model_dir="ESMDisPred-1" ;;
    2|ESMDisPred-2)   _model_dir="ESMDisPred-2" ;;
    4|ESMDisPred-DNN) _model_dir="ESMDisPred-DNN" ;;
    *)                _model_dir="ESMDisPred-2PDB" ;;  # default: model 3
  esac

  # Collect only the final ESMDisPred predictions — DisPredict3.0 also copies
  # its outputs into disorder/Dispredict3.0/, which must be excluded.
  mapfile -t caid_files < <(find "$tmp_out/disorder/$_model_dir" -name "*.caid" 2>/dev/null | sort)
  if [[ ${#caid_files[@]} -eq 0 ]]; then
    echo "ERROR: No .caid files found in disorder/$_model_dir/ for mutation size ${i}" >&2
    rm -rf "$tmp_out"
    exit 1
  fi
  cat "${caid_files[@]}" > "${OUTPUT_DIR}/${output_name}"
  rm -rf "$tmp_out"
  echo "Created ${OUTPUT_DIR}/${output_name}"
done
