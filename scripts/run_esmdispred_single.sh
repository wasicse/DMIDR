#!/usr/bin/env bash
# Run ESMDisPred on a single FASTA file via the official Docker runner.
# Must be called from the project root.
#
# Usage:
#   bash scripts/run_esmdispred_single.sh \
#     <fasta_file> <output_dir> <image> <large_models_dir> <output_stem> [model=3]
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <fasta_file> <output_dir> <image> <large_models_dir> <output_stem> [model]" >&2
  exit 1
fi

FASTA="$(realpath "$1")"
OUTPUT_DIR="$(realpath "$2")"
IMAGE_NAME="$3"
LARGE_MODELS_DIR="$(realpath "$4")"
OUTPUT_STEM="$5"
MODEL="${6:-3}"

ESMDISPRED_DIR="$(realpath "$(dirname "$0")/tools/ESMDisPred")"

if [[ ! -d "$ESMDISPRED_DIR" ]]; then
  echo "Error: ESMDisPred tool not found at $ESMDISPRED_DIR" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
tmp_out="$(mktemp -d)"

fasta_filename="$(basename "$FASTA")"
ESMpath="/opt/ESMDisPred"

# Run from inside the ESMDisPred directory so $(pwd) resolves all script mounts correctly.
pushd "$ESMDISPRED_DIR" > /dev/null

mkdir -p features

docker run --rm \
    --gpus all \
  --user "$(id -u):$(id -g)" \
  -e HOME="$ESMpath" \
  -e XDG_CACHE_HOME="$ESMpath/.cache" \
  -e TORCH_HOME="$ESMpath/largeModels" \
  -v "$FASTA":"$ESMpath/example/$fasta_filename":ro \
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
  echo "ERROR: No .caid output from ESMDisPred for $FASTA" >&2
  rm -rf "$tmp_out"
  exit 1
fi

cat "${caid_files[@]}" > "${OUTPUT_DIR}/${OUTPUT_STEM}.caid"
rm -rf "$tmp_out"
echo "Created ${OUTPUT_DIR}/${OUTPUT_STEM}.caid"
