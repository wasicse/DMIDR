from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from .io_utils import read_disorder_file


def compare_disorder_predictions(original_path: str | Path, mutant_path: str | Path) -> pd.DataFrame:
    original_df = read_disorder_file(original_path)
    mutant_df = read_disorder_file(mutant_path)

    comparison_df = pd.merge(
        original_df,
        mutant_df,
        on='Position',
        suffixes=('_orig', '_mut'),
        how='inner',
    )
    if comparison_df.empty:
        raise ValueError('No overlapping positions found between original and mutant disorder files.')

    comparison_df['DeltaDisorder'] = comparison_df['DisorderProb_mut'] - comparison_df['DisorderProb_orig']
    comparison_df['LabelChange'] = comparison_df['Disordered_mut'] - comparison_df['Disordered_orig']
    return comparison_df


def build_summary(comparison_df: pd.DataFrame) -> dict[str, float | int]:
    return {
        'total_residues': int(len(comparison_df)),
        'average_delta_disorder': float(comparison_df['DeltaDisorder'].mean()),
        'residues_increased_disorder': int((comparison_df['DeltaDisorder'] > 0).sum()),
        'residues_decreased_disorder': int((comparison_df['DeltaDisorder'] < 0).sum()),
        'label_unchanged': int((comparison_df['LabelChange'] == 0).sum()),
        'ordered_to_disordered': int((comparison_df['LabelChange'] == 1).sum()),
        'disordered_to_ordered': int((comparison_df['LabelChange'] == -1).sum()),
    }


def save_disorder_outputs(comparison_df: pd.DataFrame, output_dir: str | Path, stem: str) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    probability_csv = output_dir / f'{stem}_disorder_probability_comparison.csv'
    label_csv = output_dir / f'{stem}_disorder_label_comparison.csv'
    plot_png = output_dir / f'{stem}_disorder_score_plot.png'
    summary_txt = output_dir / f'{stem}_summary.txt'

    comparison_df[[
        'Position', 'AA_orig', 'AA_mut', 'DisorderProb_orig', 'DisorderProb_mut', 'DeltaDisorder'
    ]].to_csv(probability_csv, index=False)

    comparison_df[[
        'Position', 'AA_orig', 'AA_mut', 'Disordered_orig', 'Disordered_mut', 'LabelChange'
    ]].to_csv(label_csv, index=False)

    summary = build_summary(comparison_df)
    with summary_txt.open('w', encoding='utf-8') as handle:
        for key, value in summary.items():
            handle.write(f'{key}: {value}\n')

    plt.figure(figsize=(10, 4))
    plt.plot(comparison_df['Position'], comparison_df['DisorderProb_orig'], label='Original', marker='o')
    plt.plot(comparison_df['Position'], comparison_df['DisorderProb_mut'], label='Mutant', marker='x')
    plt.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    plt.title('Disorder Score Comparison')
    plt.xlabel('Residue Position')
    plt.ylabel('Disorder Probability')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_png, dpi=200)
    plt.close()
