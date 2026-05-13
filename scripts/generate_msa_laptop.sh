#!/usr/bin/env bash
# Run on your LOCAL LAPTOP (no firewall).
# Hits the ColabFold MSA server and saves .a3m files locally.
# Then rsync the output directory to the server for folding.
#
# Usage (from project root):
#   bash scripts/generate_msa_laptop.sh <input.fasta> [msa_output_dir]
#
# Upload results to server:
#   rsync -avz <msa_output_dir>/ <user>@<server>:<project>/data/input/msas/
#
# Then on the server run:
#   bash scripts/run_colabfold.sh data/input/msas/ results/example/alphafold
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash scripts/generate_msa_laptop.sh <input.fasta> [msa_output_dir]"
  echo ""
  echo "  input.fasta    Path to your FASTA file"
  echo "  msa_output_dir Where to save .a3m files (default: msas/)"
  exit 1
fi

FASTA="$1"
MSA_DIR="${2:-msas}"

if [[ ! -f "$FASTA" ]]; then
  echo "Error: FASTA file not found: $FASTA" >&2
  exit 1
fi

mkdir -p "$MSA_DIR"

if command -v uv &>/dev/null; then
  uv run colabfold_batch --msa-only "$FASTA" "$MSA_DIR"
elif command -v colabfold_batch &>/dev/null; then
  colabfold_batch --msa-only "$FASTA" "$MSA_DIR"
else
  echo "Error: colabfold not found. Install with: pip install 'colabfold[alphafold]'" >&2
  exit 1
fi

echo ""
echo "MSAs written to: $MSA_DIR"
echo ""
echo "Upload to server with:"
echo "  rsync -avz ${MSA_DIR}/ <user>@<server>:<project_path>/data/input/msas/"
