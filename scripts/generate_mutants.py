#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from idr_project.mutant_generator import run_mutant_generation


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Generate mutant FASTA files from a PSSM and disorder prediction.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--pssm', required=True, help='Path to ASCII PSI-BLAST PSSM file.')
    parser.add_argument('--disorder', required=True, help='Path to original .caid disorder file.')
    parser.add_argument('--output-dir', required=True, help='Directory where FASTA outputs will be written.')
    parser.add_argument('--sequence-name', required=True, help='Short name used in output filenames.')
    parser.add_argument('--mode', default='random', choices=['random', 'consecutive'],
                        help='Mutation strategy: random (scattered positions) or consecutive (sliding window blocks).')
    parser.add_argument('--min-disorder-prob', type=float, default=0.5,
                        help='Minimum disorder probability for a site to be eligible.')

    # random-mode args
    rnd = parser.add_argument_group('random mode')
    rnd.add_argument('--num-variants', type=int, default=200,
                     help='Number of random mutant variants to generate.')
    rnd.add_argument('--mutations-per-seq', type=int, default=5,
                     help='Number of mutations per variant.')
    rnd.add_argument('--min-spacing', type=int, default=5,
                     help='Minimum residue spacing between selected mutation sites.')
    rnd.add_argument('--seed', type=int, default=None,
                     help='Random seed for reproducibility (omit for non-deterministic).')

    # consecutive-mode args
    con = parser.add_argument_group('consecutive mode')
    con.add_argument('--max-mutations', type=int, default=5,
                     help='Maximum consecutive mutation block size.')
    con.add_argument('--required-ratio', type=float, default=0.8,
                     help='Minimum fraction of residues in a block that must be eligible.')

    args = parser.parse_args()

    summary = run_mutant_generation(
        pssm_path=args.pssm,
        disorder_path=args.disorder,
        output_dir=args.output_dir,
        sequence_name=args.sequence_name,
        mode=args.mode,
        num_variants=args.num_variants,
        mutations_per_seq=args.mutations_per_seq,
        min_spacing=args.min_spacing,
        seed=args.seed,
        max_mutations=args.max_mutations,
        min_disorder_prob=args.min_disorder_prob,
        required_ratio=args.required_ratio,
    )

    print(f'Mode: {args.mode}')
    print('Mutation summary:')
    for item in summary:
        print(f'  block_size={item.block_size}  mutants_created={item.mutants_created}')
    print(f'Output directory: {Path(args.output_dir).resolve()}')


if __name__ == '__main__':
    main()
