#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description='Overlay two time-series analysis CSV files.')
    parser.add_argument('analysis_a')
    parser.add_argument('analysis_b')
    parser.add_argument('outdir')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dfa = pd.read_csv(args.analysis_a)
    dfb = pd.read_csv(args.analysis_b)

    plt.figure()
    plt.plot(dfa['time_ps'], dfa['Rg_nm'], label='Run A')
    plt.plot(dfb['time_ps'], dfb['Rg_nm'], label='Run B')
    plt.xlabel('Time (ps)')
    plt.ylabel('Rg (nm)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / 'Rg_overlay.png', dpi=180)
    plt.close()

    plt.figure()
    plt.plot(dfa['time_ps'], dfa['end_to_end_nm'], label='Run A')
    plt.plot(dfb['time_ps'], dfb['end_to_end_nm'], label='Run B')
    plt.xlabel('Time (ps)')
    plt.ylabel('End-to-end distance (nm)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / 'EndToEnd_overlay.png', dpi=180)
    plt.close()

    summary = pd.DataFrame({
        'metric': ['Rg_nm', 'end_to_end_nm'],
        'A_mean': [dfa['Rg_nm'].mean(), dfa['end_to_end_nm'].mean()],
        'B_mean': [dfb['Rg_nm'].mean(), dfb['end_to_end_nm'].mean()],
        'A_std': [dfa['Rg_nm'].std(), dfa['end_to_end_nm'].std()],
        'B_std': [dfb['Rg_nm'].std(), dfb['end_to_end_nm'].std()],
    })
    summary.to_csv(outdir / 'summary_AB.csv', index=False)
    print(f'Wrote outputs to {outdir}')


if __name__ == '__main__':
    main()
