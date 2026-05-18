#!/usr/bin/env bash
# Waits for the running 1res container to finish, collects its output,
# then re-runs run_all.sh steps 4-5 for the remaining mutation sizes.
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
CONTAINER_ID="a6049c881f6a"
TMP_OUT="/tmp/tmp.lHs76CGN0c"
OUTPUT_DIR="$PROJECT_ROOT/outputs/Mur18B-PA/dispred"
CAID_OUT="$OUTPUT_DIR/Mur18B-PA_mutants_1res.caid"

echo "[$(date '+%H:%M:%S')] Waiting for container $CONTAINER_ID to finish..."
docker wait "$CONTAINER_ID" || true

echo "[$(date '+%H:%M:%S')] Container done. Collecting 1res output..."
mkdir -p "$OUTPUT_DIR"
mapfile -t caid_files < <(find "$TMP_OUT" -name "*.caid" | sort)

if [[ ${#caid_files[@]} -eq 0 ]]; then
    echo "ERROR: No .caid files found in $TMP_OUT" >&2
    ls -la "$TMP_OUT" >&2
    exit 1
fi

cat "${caid_files[@]}" > "$CAID_OUT"
echo "[$(date '+%H:%M:%S')] Saved → $CAID_OUT"
rm -rf "$TMP_OUT"

echo "[$(date '+%H:%M:%S')] Continuing pipeline (steps 4-5 for remaining mutation sizes)..."
echo "4 5" | bash "$PROJECT_ROOT/run_all.sh" "$PROJECT_ROOT/configs/local.env"
