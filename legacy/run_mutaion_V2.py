import pandas as pd
import re
import random
from io import StringIO

# === PSSM Reader ===
def read_blast_pssm(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    pssm_lines = []
    pssm_start_found = False

    for line in lines:
        if re.match(r"\s*\d+\s+[A-Z]\s+-?\d", line):
            pssm_start_found = True
            pssm_lines.append(line.strip())
        elif pssm_start_found:
            if line.strip() == "":
                break
            pssm_lines.append(line.strip())

    if not pssm_lines:
        raise ValueError("No PSSM data found in the file.")

    cols = ["Position", "Residue"] + [
        "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    ] + [f"Perc_{aa}" for aa in [
        "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    ]] + ["Information", "RelativeWeight"]

    data_str = "\n".join(pssm_lines)
    df = pd.read_csv(StringIO(data_str), delim_whitespace=True, names=cols)

    return df


# === Disorder File Reader ===
def read_disorder_file(file_path):
    disorder_data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>') or line.strip() == '':
                continue
            parts = line.strip().split()
            if len(parts) == 4:
                seq_no = int(parts[0])
                aa = parts[1]
                prob = float(parts[2])
                flag = int(parts[3])
                disorder_data.append([seq_no, aa, prob, flag])

    df = pd.DataFrame(disorder_data, columns=['SeqNo', 'AA', 'DisorderProb', 'Disordered'])
    return df


# === Multi-Mutation Generator ===
import os

def generate_multiple_mutants(pssm_df, disorder_df, num_variants=5, mutations_per_seq=5, output_fasta="mutated_sequences.fasta"):
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    mutation_log_path = os.path.join(os.path.dirname(output_fasta), "mutation_log.txt")

    # Clear mutation log at the beginning
    open(mutation_log_path, "w").close()

    if 'SeqNo' in disorder_df.columns:
        disorder_df = disorder_df.rename(columns={'SeqNo': 'Position'})

    merged_df = pd.merge(pssm_df, disorder_df, on='Position')
    disorder_promoting = set(['E', 'K', 'R', 'Q', 'S', 'P', 'G'])
    order_promoting = ['L', 'I', 'F', 'V', 'Y', 'W', 'M', 'C']

    original_sequence = list(merged_df['Residue'])

    with open(output_fasta, 'w') as fasta_file:
        for variant_idx in range(1, num_variants + 1):
            candidates = merged_df[
                (merged_df['Disordered'] == 1) &
                (merged_df['Residue'].isin(disorder_promoting)) &
                (merged_df['DisorderProb'] > 0.5)
            ].copy()

            candidates = candidates.sort_values(by=['DisorderProb', 'Information'], ascending=[False, True])
            top_k = candidates.head(100)
            positions = top_k['Position'].tolist()
            random.shuffle(positions)

            selected_positions = []
            for pos in positions:
                if all(abs(pos - sel) >= 5 for sel in selected_positions):
                    selected_positions.append(pos)
                if len(selected_positions) >= mutations_per_seq * 3:
                    break

            mutated_sequence = original_sequence.copy()
            actual_mutations = 0
            mutation_details = []

            for pos in selected_positions:
                if actual_mutations >= mutations_per_seq:
                    break

                row = merged_df[merged_df['Position'] == pos].iloc[0]
                original_aa = row['Residue']
                scores = row[order_promoting].dropna().astype(float).sort_values(ascending=False)

                for alt in scores.index:
                    if alt != original_aa:
                        mutated_sequence[pos - 1] = alt
                        mutation_details.append(f"{original_aa}{pos}{alt}")
                        actual_mutations += 1
                        break

            if actual_mutations < mutations_per_seq:
                print(f"⚠️ Warning: Only {actual_mutations} mutations made for Mutant_{variant_idx} (expected {mutations_per_seq})")

            header = f">Mutant_{variant_idx}"
            with open(mutation_log_path, "a") as log_file:
                print(f"{header} | Mutations: {';'.join(mutation_details)}", file=log_file)

            fasta_file.write(header + "\n")
            for i in range(0, len(mutated_sequence), 60):
                fasta_file.write(''.join(mutated_sequence[i:i + 60]) + "\n")

    print(f"\n✅ Generated {num_variants} mutated sequences (each with exactly {mutations_per_seq} changes if possible) in: {output_fasta}")



# === Main Runner ===
if __name__ == "__main__":
    # === File Paths ===
    pssm_path = 'fonPB_O.pssm'
    disorder_path = './Dispred/fonPB_O.dispred'
    output_fasta_path = './Fasta/fonPB_MultiMutants.fasta'

    # === Load Data ===
    pssm_df = read_blast_pssm(pssm_path)
    disorder_df = read_disorder_file(disorder_path)

    # === Generate Multiple Mutants ===
    generate_multiple_mutants(pssm_df, disorder_df, num_variants=200, mutations_per_seq=5, output_fasta=output_fasta_path)
