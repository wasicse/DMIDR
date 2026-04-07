#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from idr_project.mutant_generator import run_mutant_generation


def main() -> None:
    parser = argparse.ArgumentParser(description='Generate mutant FASTA files from a PSSM and Dispredict result.')
    parser.add_argument('--pssm', required=True, help='Path to ASCII PSI-BLAST PSSM file.')
    parser.add_argument('--disorder', required=True, help='Path to original .dispred file.')
    parser.add_argument('--output-dir', required=True, help='Directory where FASTA outputs will be written.')
    parser.add_argument('--sequence-name', required=True, help='Short name used in output filenames.')
    parser.add_argument('--max-mutations', type=int, default=5, help='Maximum consecutive mutation block size.')
    parser.add_argument('--min-disorder-prob', type=float, default=0.5, help='Minimum disorder probability to consider a site eligible.')
    parser.add_argument('--required-ratio', type=float, default=0.8, help='Minimum fraction of residues in a block that must be eligible.')
    args = parser.parse_args()

    summary = run_mutant_generation(
        pssm_path=args.pssm,
        disorder_path=args.disorder,
        output_dir=args.output_dir,
        sequence_name=args.sequence_name,
        max_mutations=args.max_mutations,
        min_disorder_prob=args.min_disorder_prob,
        required_ratio=args.required_ratio,
    )

    print('Mutation summary:')
    for item in summary:
        print(f'  block_size={item.block_size} mutants_created={item.mutants_created}')
    print(f'Output directory: {Path(args.output_dir).resolve()}')


if __name__ == '__main__':
    main()
