from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .io_utils import read_blast_pssm, read_disorder_file, write_fasta


DISORDER_PROMOTING = {'E', 'K', 'R', 'Q', 'S', 'P', 'G'}
ORDER_PROMOTING = ['L', 'I', 'F', 'V', 'Y', 'W', 'M', 'C']


@dataclass
class MutationSummary:
    block_size: int
    mutants_created: int


def merge_pssm_and_disorder(pssm_df: pd.DataFrame, disorder_df: pd.DataFrame) -> pd.DataFrame:
    disorder_df = disorder_df.rename(columns={'AA': 'DisorderAA'})
    merged = pd.merge(pssm_df, disorder_df, on='Position', how='inner')
    if merged.empty:
        raise ValueError('PSSM and disorder tables could not be merged on Position.')
    return merged


def generate_consecutive_mutants(
    pssm_df: pd.DataFrame,
    disorder_df: pd.DataFrame,
    max_mutations: int = 5,
    min_disorder_prob: float = 0.5,
    required_ratio: float = 0.8,
) -> dict[int, list[tuple[str, str]]]:
    """Generate mutant FASTA records grouped by mutation block size.

    Returns a dict like {1: [(header, sequence), ...], 2: [...]}.
    """
    merged_df = merge_pssm_and_disorder(pssm_df, disorder_df)
    original_sequence = list(merged_df['Residue'])
    sequence_length = len(original_sequence)

    mutants_by_size: dict[int, list[tuple[str, str]]] = defaultdict(list)
    mutant_counter: dict[int, int] = defaultdict(int)

    for block_size in range(1, max_mutations + 1):
        for start_pos in range(1, sequence_length - block_size + 2):
            block_positions = list(range(start_pos, start_pos + block_size))
            rows = merged_df[merged_df['Position'].isin(block_positions)]
            if len(rows) < block_size:
                continue

            valid_count = sum(
                (row['Disordered'] == 1 and row['Residue'] in DISORDER_PROMOTING and row['DisorderProb'] > min_disorder_prob)
                for _, row in rows.iterrows()
            )
            if valid_count / block_size < required_ratio:
                continue

            mutated_sequence = original_sequence.copy()
            mutation_details: list[str] = []
            success = True

            for _, row in rows.iterrows():
                pos = int(row['Position'])
                original_aa = row['Residue']
                candidate_scores = row[ORDER_PROMOTING].dropna().astype(float).sort_values(ascending=False)
                alt_aa = next((aa for aa in candidate_scores.index if aa != original_aa), None)
                if alt_aa is None:
                    success = False
                    break
                mutated_sequence[pos - 1] = alt_aa
                mutation_details.append(f'{original_aa}{pos}{alt_aa}')

            if not success:
                continue

            mutant_counter[block_size] += 1
            header = f"Mutant_{mutant_counter[block_size]}_Res{block_size}_" + '_'.join(mutation_details)
            mutants_by_size[block_size].append((header, ''.join(mutated_sequence)))

    return dict(mutants_by_size)


def save_mutants(mutants_by_size: dict[int, list[tuple[str, str]]], output_dir: str | Path, prefix: str) -> list[MutationSummary]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary: list[MutationSummary] = []
    log_lines: list[str] = []

    for block_size in sorted(mutants_by_size):
        records = mutants_by_size[block_size]
        fasta_path = output_dir / f'{prefix}_mutants_{block_size}res.fasta'
        write_fasta(fasta_path, records)
        log_lines.extend([f'>{header}' for header, _ in records])
        summary.append(MutationSummary(block_size=block_size, mutants_created=len(records)))

    (output_dir / f'{prefix}_mutation_log.txt').write_text('\n'.join(log_lines) + ('\n' if log_lines else ''), encoding='utf-8')
    return summary


def run_mutant_generation(
    pssm_path: str | Path,
    disorder_path: str | Path,
    output_dir: str | Path,
    sequence_name: str,
    max_mutations: int = 5,
    min_disorder_prob: float = 0.5,
    required_ratio: float = 0.8,
) -> list[MutationSummary]:
    pssm_df = read_blast_pssm(pssm_path)
    disorder_df = read_disorder_file(disorder_path)
    mutants = generate_consecutive_mutants(
        pssm_df=pssm_df,
        disorder_df=disorder_df,
        max_mutations=max_mutations,
        min_disorder_prob=min_disorder_prob,
        required_ratio=required_ratio,
    )
    return save_mutants(mutants, output_dir=output_dir, prefix=sequence_name)
