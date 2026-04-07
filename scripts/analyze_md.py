#!/usr/bin/env python3
from __future__ import annotations

import argparse

from idr_project.md_analysis import run_md_analysis


def main() -> None:
    parser = argparse.ArgumentParser(description='Analyze WT and mutant MD trajectories.')
    parser.add_argument('--wt-top', required=True)
    parser.add_argument('--wt-traj', required=True)
    parser.add_argument('--mut-top', required=True)
    parser.add_argument('--mut-traj', required=True)
    parser.add_argument('--disorder-table', required=True, help='CSV/TSV with columns resid,disorder_prob')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--cutoff-nm', type=float, default=0.8)
    parser.add_argument('--exclude-neighbors', type=int, default=1)
    args = parser.parse_args()

    run_md_analysis(
        wt_top=args.wt_top,
        wt_traj=args.wt_traj,
        mut_top=args.mut_top,
        mut_traj=args.mut_traj,
        disorder_path=args.disorder_table,
        outdir=args.outdir,
        cutoff_nm=args.cutoff_nm,
        exclude_neighbors=args.exclude_neighbors,
    )
    print(f'Wrote MD analysis results to {args.outdir}')


if __name__ == '__main__':
    main()
