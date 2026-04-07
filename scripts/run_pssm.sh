#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: bash scripts/run_pssm.sh <query_fasta> <blast_db> <output_pssm> [threads] [iterations]"
  exit 1
fi

QUERY_FASTA="$1"
BLAST_DB="$2"
OUTPUT_PSSM="$3"
THREADS="${4:-8}"
ITERATIONS="${5:-3}"

psiblast \
  -query "$QUERY_FASTA" \
  -db "$BLAST_DB" \
  -num_iterations "$ITERATIONS" \
  -out_ascii_pssm "$OUTPUT_PSSM" \
  -num_threads "$THREADS"
