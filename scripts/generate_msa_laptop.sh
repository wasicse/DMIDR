#!/usr/bin/env bash
# Run on your LOCAL LAPTOP (no firewall).
# Hits the ColabFold MSA server and saves .a3m files locally.
# Then rsync the output directory to the server for folding.
#
# Usage:
#   bash scripts/generate_msa_laptop.sh <input.fasta> [msa_output_dir]
#
# Upload results to server:
#   rsync -avz <msa_output_dir>/ <user>@<server>:<project>/data/input/msas/
#
# Then on the server run:
#   bash scripts/run_colabfold.sh data/input/msas/ results/example/alphafold
set -euo pipefail

FASTA="${1:?Usage: $0 <input.fasta> [msa_output_dir]}"
MSA_DIR="${2:-msas}"

if [[ ! -f "$FASTA" ]]; then
  echo "Error: FASTA file not found: $FASTA" >&2
  exit 1
fi

mkdir -p "$MSA_DIR"

colabfold_batch --msa-only "$FASTA" "$MSA_DIR"

echo "MSAs written to: $MSA_DIR"
echo ""
echo "Upload to server with:"
echo "  rsync -avz ${MSA_DIR}/ <user>@<server>:<project_path>/data/input/msas/"
