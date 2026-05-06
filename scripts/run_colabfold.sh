#!/usr/bin/env bash
set -euo pipefail

FASTA_DEFAULT="/home/mkabir3/Research/47_Dr_Raj_Colab/DMIDR/input/example.fasta"
OUTDIR_DEFAULT="result/alphafold"

fasta_path="${1:-$FASTA_DEFAULT}"
output_dir="${2:-$OUTDIR_DEFAULT}"
shift $(( $# > 0 ? 1 : 0 )) || true
shift $(( $# > 0 ? 1 : 0 )) || true

if ! command -v uv >/dev/null 2>&1; then
  echo "Error: uv not found. Install it from https://docs.astral.sh/uv/" >&2
  exit 1
fi

# Accept either a FASTA file or a directory of pre-computed .a3m MSAs as input.
if [[ -d "$fasta_path" ]]; then
  echo "Input is a directory — using pre-computed MSAs from: $fasta_path"
elif [[ ! -f "$fasta_path" ]]; then
  echo "Error: input not found (expected FASTA file or MSA directory): $fasta_path" >&2
  exit 1
fi

mkdir -p "$output_dir"

# Pass any extra args after the first two directly to colabfold_batch.
uv run colabfold_batch "$@" "$fasta_path" "$output_dir"

echo "ColabFold run complete. Outputs in: $output_dir"
