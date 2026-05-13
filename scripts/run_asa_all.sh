#!/usr/bin/env bash
# Compute per-residue ASA for every PDB file found under a directory.
# No external binaries required — uses mdtraj via the uv-managed environment.
#
# Usage:
#   bash scripts/run_asa_all.sh <pdb_dir> <results_dir> [extra calculate_asa.py flags]
#
# Example:
#   bash scripts/run_asa_all.sh results/example/alphafold results/example/asa
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

PDB_DIR="${1:?Usage: $0 <pdb_dir> <results_dir> [extra flags]}"
RESULTS_DIR="${2:?Usage: $0 <pdb_dir> <results_dir> [extra flags]}"
shift 2

if ! command -v uv &>/dev/null; then
  echo "Error: uv not found. Install from https://docs.astral.sh/uv/" >&2
  exit 1
fi

mapfile -t PDBS < <(find "$PDB_DIR" -name "*.pdb" | sort)

if [[ ${#PDBS[@]} -eq 0 ]]; then
  echo "No PDB files found in: $PDB_DIR" >&2
  exit 1
fi

echo "Found ${#PDBS[@]} PDB file(s). Writing CSVs to: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

for pdb in "${PDBS[@]}"; do
  stem="$(basename "$pdb" .pdb)"
  out="$RESULTS_DIR/${stem}_asa.csv"
  echo "  Processing: $stem"
  uv run --project "$PROJECT_ROOT" python "$PROJECT_ROOT/scripts/calculate_asa.py" \
    --pdb "$pdb" \
    --out "$out" \
    "$@"
done

echo ""
echo "Done. CSVs written to: $RESULTS_DIR"
