#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: bash scripts/run_dispredict.sh <fasta_dir> <output_dir> <image_name> <container_prefix> [max_mutations]"
  exit 1
fi

FASTA_DIR="$1"
OUTPUT_DIR="$2"
IMAGE_NAME="$3"
CONTAINER_PREFIX="$4"
MAX_MUTATIONS="${5:-5}"

mkdir -p "$OUTPUT_DIR"

for i in $(seq 1 "$MAX_MUTATIONS"); do
  file=$(find "$FASTA_DIR" -maxdepth 1 -type f -name "*_mutants_${i}res.fasta" | head -n 1)
  if [[ -z "$file" ]]; then
    echo "Skipping mutation size ${i}: no FASTA file found"
    continue
  fi

  container_name="${CONTAINER_PREFIX}_${i}"
  output_name="$(basename "${file%.fasta}").dispred"

  docker stop "$container_name" >/dev/null 2>&1 || true
  docker rm "$container_name" >/dev/null 2>&1 || true
  docker run -ti -d --name "$container_name" "$IMAGE_NAME" >/dev/null

  docker cp "$file" "$container_name":/opt/Dispredict3.0/example/sample.fasta
  docker exec -i "$container_name" /bin/bash -lc \
    "source /opt/Dispredict3.0/.venv/bin/activate && \
     /opt/Dispredict3.0/.venv/bin/python /opt/Dispredict3.0/script/Dispredict3.0.py \
       -f /opt/Dispredict3.0/example/sample.fasta \
       -o /opt/Dispredict3.0/output/"
  docker cp "$container_name":/opt/Dispredict3.0/output/sample_disPred.txt "$OUTPUT_DIR/$output_name"

  echo "Created $OUTPUT_DIR/$output_name"
done
