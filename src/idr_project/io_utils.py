from __future__ import annotations

import re
from io import StringIO
from pathlib import Path

import pandas as pd


PSSM_SCORE_COLUMNS = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
]

PSSM_PERCENT_COLUMNS = [f'Perc_{aa}' for aa in PSSM_SCORE_COLUMNS]

DISORDER_THRESHOLD = 0.5


def read_blast_pssm(file_path: str | Path) -> pd.DataFrame:
    """Read an ASCII PSI-BLAST PSSM file into a dataframe."""
    file_path = Path(file_path)
    lines = file_path.read_text(encoding='utf-8', errors='ignore').splitlines()

    pssm_lines: list[str] = []
    pssm_started = False
    for line in lines:
        if re.match(r"\s*\d+\s+[A-Z]\s+-?\d", line):
            pssm_started = True
            pssm_lines.append(line.strip())
        elif pssm_started:
            if not line.strip():
                break
            pssm_lines.append(line.strip())

    if not pssm_lines:
        raise ValueError(f'No PSSM data found in {file_path}')

    columns = ['Position', 'Residue'] + PSSM_SCORE_COLUMNS + PSSM_PERCENT_COLUMNS + ['Information', 'RelativeWeight']
    df = pd.read_csv(StringIO('\n'.join(pssm_lines)), sep=r'\s+', names=columns, engine='python')
    df['Position'] = df['Position'].astype(int)
    return df


def read_disorder_file(file_path: str | Path) -> pd.DataFrame:
    """Read an ESMDisPred .caid or legacy Dispredict .dispred file into a dataframe.

    .caid format (3 columns, no header lines):
        position  amino_acid  disorder_probability

    .dispred format (4 columns, with >header lines):
        position  amino_acid  disorder_probability  disordered_flag

    The binary Disordered column is read directly for .dispred files and
    computed as prob >= 0.5 for .caid files.
    """
    file_path = Path(file_path)
    disorder_data: list[list[object]] = []

    with file_path.open('r', encoding='utf-8', errors='ignore') as handle:
        seen_header = False
        for line in handle:
            if line.startswith('>'):
                if seen_header:
                    break  # stop at second sequence — only use first chain
                seen_header = True
                continue
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) == 3:
                prob = float(parts[2])
                disorder_data.append([int(parts[0]), parts[1], prob, int(prob >= DISORDER_THRESHOLD)])
            elif len(parts) == 4:
                disorder_data.append([int(parts[0]), parts[1], float(parts[2]), int(parts[3])])

    if not disorder_data:
        raise ValueError(f'No disorder rows found in {file_path}')

    return pd.DataFrame(disorder_data, columns=['Position', 'AA', 'DisorderProb', 'Disordered'])


def write_fasta(path: str | Path, records: list[tuple[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8') as handle:
        for header, sequence in records:
            handle.write(f'>{header}\n')
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + '\n')


def load_disorder_probability_table(path: str | Path) -> pd.DataFrame:
    """Load a CSV/TSV file with columns resid, disorder_prob.

    Also accepts the pipeline's disorder probability comparison CSV which uses
    Position / DisorderProb_orig column names.
    """
    path = Path(path)
    with path.open('r', encoding='utf-8', errors='ignore') as handle:
        first_line = handle.readline()
    sep = ',' if ',' in first_line else '\t' if '\t' in first_line else r'\s+'
    df = pd.read_csv(path, sep=sep, engine='python')
    # Normalise alternative column names from the comparison CSV
    df = df.rename(columns={'Position': 'resid', 'DisorderProb_orig': 'disorder_prob'})
    required = {'resid', 'disorder_prob'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'{path} is missing columns: {sorted(missing)}')
    df = df.copy()
    df['resid'] = df['resid'].astype(int)
    df['disorder_prob'] = df['disorder_prob'].astype(float)
    return df.sort_values('resid')
