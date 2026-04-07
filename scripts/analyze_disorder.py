#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from idr_project.disorder_analysis import compare_disorder_predictions, save_disorder_outputs, build_summary


def main() -> None:
    parser = argparse.ArgumentParser(description='Compare original and mutant Dispredict outputs.')
    parser.add_argument('--original', required=True, help='Original sequence .dispred file')
    parser.add_argument('--mutant', required=True, help='Mutant sequence .dispred file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--label', default='comparison', help='Prefix used for output filenames')
    args = parser.parse_args()

    comparison = compare_disorder_predictions(args.original, args.mutant)
    save_disorder_outputs(comparison, args.output_dir, args.label)
    summary = build_summary(comparison)

    print('Disorder comparison summary:')
    for key, value in summary.items():
        print(f'  {key}: {value}')
    print(f'Output directory: {Path(args.output_dir).resolve()}')


if __name__ == '__main__':
    main()
