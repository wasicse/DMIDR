#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG_FILE="${1:-$PROJECT_ROOT/configs/local.env}"

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  echo "Create it from: $PROJECT_ROOT/configs/example.env" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"

: "${SEQ_NAME:?SEQ_NAME is required in the config file}"
: "${INPUT_FASTA:?INPUT_FASTA is required in the config file}"
: "${BLAST_DB:?BLAST_DB is required in the config file}"
: "${DISPREDICT_IMAGE:?DISPREDICT_IMAGE is required in the config file}"
MAX_MUTATIONS="${MAX_MUTATIONS:-5}"
CONTAINER_PREFIX="${CONTAINER_PREFIX:-$SEQ_NAME}"
SKIP_PSSM="${SKIP_PSSM:-0}"
SKIP_ORIGINAL_DISPREDICT="${SKIP_ORIGINAL_DISPREDICT:-0}"
SKIP_MUTANT_GENERATION="${SKIP_MUTANT_GENERATION:-0}"
SKIP_MUTANT_DISPREDICT="${SKIP_MUTANT_DISPREDICT:-0}"

PSSM_OUT="${PSSM_OUT:-$PROJECT_ROOT/data/intermediate/pssm/${SEQ_NAME}.pssm}"
DISPRED_DIR="${DISPRED_DIR:-$PROJECT_ROOT/data/intermediate/dispred}"
FASTA_OUT="${FASTA_OUT:-$PROJECT_ROOT/data/intermediate/fasta/${SEQ_NAME}}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_ROOT/results/${SEQ_NAME}}"
ORIGINAL_DISPRED="${ORIGINAL_DISPRED:-$DISPRED_DIR/${SEQ_NAME}_original.dispred}"

mkdir -p "$(dirname "$PSSM_OUT")" "$DISPRED_DIR" "$FASTA_OUT" "$RESULTS_DIR"

step() {
  echo
  echo "============================================================"
  echo "$1"
  echo "============================================================"
}

run_original_dispredict() {
  local container_name="${CONTAINER_PREFIX}_orig"
  docker stop "$container_name" >/dev/null 2>&1 || true
  docker rm "$container_name" >/dev/null 2>&1 || true
  docker run -ti -d --name "$container_name" "$DISPREDICT_IMAGE" >/dev/null
  docker cp "$INPUT_FASTA" "$container_name":/opt/Dispredict3.0/example/sample.fasta
  docker exec -i "$container_name" /bin/bash -lc \
    "source /opt/Dispredict3.0/.venv/bin/activate && \
     /opt/Dispredict3.0/.venv/bin/python /opt/Dispredict3.0/script/Dispredict3.0.py \
     -f /opt/Dispredict3.0/example/sample.fasta -o /opt/Dispredict3.0/output/"
  docker cp "$container_name":/opt/Dispredict3.0/output/sample_disPred.txt "$ORIGINAL_DISPRED"
  docker stop "$container_name" >/dev/null 2>&1 || true
  docker rm "$container_name" >/dev/null 2>&1 || true
}

if [[ "$SKIP_PSSM" != "1" ]]; then
  step "[1/4] Generate PSSM"
  bash "$PROJECT_ROOT/scripts/run_pssm.sh" "$INPUT_FASTA" "$BLAST_DB" "$PSSM_OUT"
else
  step "[1/4] Skipping PSSM generation"
fi

if [[ "$SKIP_ORIGINAL_DISPREDICT" != "1" ]]; then
  step "[2/4] Run Dispredict for original sequence"
  run_original_dispredict
else
  step "[2/4] Skipping original Dispredict run"
fi

if [[ "$SKIP_MUTANT_GENERATION" != "1" ]]; then
  step "[3/4] Generate mutant FASTA files"
  PYTHONPATH="$PROJECT_ROOT/src" python "$PROJECT_ROOT/scripts/generate_mutants.py" \
    --pssm "$PSSM_OUT" \
    --disorder "$ORIGINAL_DISPRED" \
    --output-dir "$FASTA_OUT" \
    --sequence-name "$SEQ_NAME" \
    --max-mutations "$MAX_MUTATIONS"
else
  step "[3/4] Skipping mutant FASTA generation"
fi

if [[ "$SKIP_MUTANT_DISPREDICT" != "1" ]]; then
  step "[4/4] Run Dispredict for mutant FASTA files"
  bash "$PROJECT_ROOT/scripts/run_dispredict.sh" \
    "$FASTA_OUT" \
    "$DISPRED_DIR" \
    "$DISPREDICT_IMAGE" \
    "$CONTAINER_PREFIX" \
    "$MAX_MUTATIONS"
else
  step "[4/4] Skipping mutant Dispredict runs"
fi

echo
printf 'Pipeline finished for %s\n' "$SEQ_NAME"
printf 'Original Dispredict: %s\n' "$ORIGINAL_DISPRED"
printf 'Mutant FASTA dir:    %s\n' "$FASTA_OUT"
printf 'Dispredict dir:      %s\n' "$DISPRED_DIR"
printf 'Results dir:         %s\n' "$RESULTS_DIR"
echo
cat <<MSG
Next examples:
  PYTHONPATH=src python scripts/analyze_disorder.py \
    --original "$ORIGINAL_DISPRED" \
    --mutant "$DISPRED_DIR/${SEQ_NAME}_mutants_3res.dispred" \
    --output-dir "$RESULTS_DIR/disorder_3res" \
    --label "${SEQ_NAME}_3res"
MSG
