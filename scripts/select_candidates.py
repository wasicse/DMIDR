#!/usr/bin/env python3
"""Select top-N mutant candidates ranked by average disorder reduction vs wild-type.

Ranking metrics (matching paper methodology):
  1. avg_disorder_reduction  — mean(DisorderProb_orig - DisorderProb_mut) across all residues
  2. disordered_to_ordered   — residues where binary label flips 1 (disordered) → 0 (ordered)

Reads all mutants_Nres.caid files (already computed by ESMDisPred in step 4),
scores every mutant, and writes the top-N as individual FASTA files.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_caid_multi(path: Path) -> dict[str, list[tuple[int, str, float, int]]]:
    """Parse a multi-sequence .caid file.

    Returns {header: [(pos, aa, prob, label), ...]}
    where label is 1=disordered, 0=ordered.
    """
    sequences: dict[str, list[tuple[int, str, float, int]]] = {}
    current_header: str | None = None
    current_residues: list[tuple[int, str, float, int]] = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = current_residues
                current_header = line[1:]
                current_residues = []
            else:
                parts = line.split()
                if len(parts) >= 3:
                    pos  = int(parts[0])
                    aa   = parts[1]
                    prob = float(parts[2])
                    label = int(parts[3]) if len(parts) >= 4 else int(prob > 0.5)
                    current_residues.append((pos, aa, prob, label))

    if current_header is not None:
        sequences[current_header] = current_residues

    return sequences


def parse_caid_single(path: Path) -> list[tuple[int, str, float, int]]:
    """Parse a single-sequence .caid file into [(pos, aa, prob, label), ...]."""
    residues: list[tuple[int, str, float, int]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                pos   = int(parts[0])
                aa    = parts[1]
                prob  = float(parts[2])
                label = int(parts[3]) if len(parts) >= 4 else int(prob > 0.5)
                residues.append((pos, aa, prob, label))
    return residues


def parse_fasta(path: Path) -> dict[str, str]:
    """Parse a FASTA file into {header: sequence}."""
    sequences: dict[str, str] = {}
    current_header: str | None = None
    current_seq: list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)

    return sequences


def score_mutants(
    orig_residues: list[tuple[int, str, float, int]],
    dispred_dir: Path,
    sequence_name: str,
    max_block_size: int,
) -> list[dict]:
    """Compute disorder reduction scores for every mutant across all block sizes."""
    orig_prob  = {pos: prob  for pos, aa, prob, label in orig_residues}
    orig_label = {pos: label for pos, aa, prob, label in orig_residues}
    seq_length = len(orig_residues)

    rows = []
    for block_size in range(1, max_block_size + 1):
        caid_path = dispred_dir / f'{sequence_name}_mutants_{block_size}res.caid'
        if not caid_path.exists():
            continue

        for header, residues in parse_caid_multi(caid_path).items():
            if len(residues) != seq_length:
                continue

            total_reduction = sum(
                orig_prob.get(pos, 0.0) - prob
                for pos, aa, prob, label in residues
            )
            disordered_to_ordered = sum(
                1 for pos, aa, prob, label in residues
                if orig_label.get(pos, 0) == 1 and label == 0
            )
            ordered_to_disordered = sum(
                1 for pos, aa, prob, label in residues
                if orig_label.get(pos, 0) == 0 and label == 1
            )

            rows.append({
                'mutant_id':             header,
                'block_size':            block_size,
                'avg_disorder_reduction': total_reduction / seq_length,
                'total_disorder_reduction': total_reduction,
                'disordered_to_ordered': disordered_to_ordered,
                'ordered_to_disordered': ordered_to_disordered,
            })

    return rows


def write_candidate_fastas(
    ranked_rows: list[dict],
    mutants_dir: Path,
    output_dir: Path,
    top_n: int,
) -> list[Path]:
    """Write top-N candidate sequences as individual FASTA files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_lookup: dict[str, str] = {}
    for fasta_path in sorted(mutants_dir.glob('*.fasta')):
        fasta_lookup.update(parse_fasta(fasta_path))

    written: list[Path] = []
    for rank, row in enumerate(ranked_rows[:top_n], start=1):
        header   = row['mutant_id']
        sequence = fasta_lookup.get(header)
        if sequence is None:
            print(f'  [WARN] Sequence not found in FASTA for: {header}')
            continue

        safe_name = header.replace(' ', '_').replace('/', '_')
        out_path  = output_dir / f'candidate_{rank}_{safe_name}.fasta'
        with open(out_path, 'w') as fh:
            fh.write(f'>{header}\n')
            for i in range(0, len(sequence), 60):
                fh.write(sequence[i:i + 60] + '\n')

        written.append(out_path)
        print(f'  [{rank}] {header}')
        print(f'       avg_disorder_reduction = {row["avg_disorder_reduction"]:.4f}')
        print(f'       disordered_to_ordered  = {row["disordered_to_ordered"]}')
        print(f'       → {out_path.name}')

    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Select top-N mutants by average disorder reduction vs wild-type.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--dispred-dir',    required=True,
                        help='Directory containing original and mutant .caid files.')
    parser.add_argument('--mutants-dir',    required=True,
                        help='Directory containing mutant FASTA files.')
    parser.add_argument('--output-dir',     required=True,
                        help='Directory to write candidate FASTA files.')
    parser.add_argument('--sequence-name',  required=True,
                        help='Sequence name prefix used in file names.')
    parser.add_argument('--top-n',          type=int, default=3,
                        help='Number of top candidates to select.')
    parser.add_argument('--max-block-size', type=int, default=5,
                        help='Maximum block size to consider.')
    args = parser.parse_args()

    dispred_dir = Path(args.dispred_dir)
    mutants_dir = Path(args.mutants_dir)
    output_dir  = Path(args.output_dir)

    orig_caid = dispred_dir / f'{args.sequence_name}_original.caid'
    if not orig_caid.exists():
        print(f'ERROR: original caid not found: {orig_caid}')
        raise SystemExit(1)

    print(f'Scoring mutants for: {args.sequence_name}')
    orig_residues = parse_caid_single(orig_caid)

    rows = score_mutants(orig_residues, dispred_dir, args.sequence_name, args.max_block_size)
    if not rows:
        print('ERROR: No mutant caid files found.')
        raise SystemExit(1)

    df = (
        pd.DataFrame(rows)
        .sort_values('avg_disorder_reduction', ascending=False)
        .reset_index(drop=True)
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / 'candidate_scores.csv'
    df.to_csv(summary_path, index=False)
    print(f'Scored {len(df)} mutants — full ranking saved to {summary_path}')
    print(f'\nTop {args.top_n} candidates (ranked by avg disorder reduction):')
    print(f'{"Rank":<5} {"MutantID":<55} {"Avg Reduction":>14} {"Dis→Ord":>8}')
    print('-' * 85)
    for rank, row in enumerate(df.head(args.top_n).to_dict('records'), 1):
        print(f'{rank:<5} {row["mutant_id"]:<55} {row["avg_disorder_reduction"]:>14.4f} {row["disordered_to_ordered"]:>8}')

    print()
    written = write_candidate_fastas(df.to_dict('records'), mutants_dir, output_dir, args.top_n)
    print(f'\nWrote {len(written)} candidate FASTA(s) to {output_dir}')


if __name__ == '__main__':
    main()
