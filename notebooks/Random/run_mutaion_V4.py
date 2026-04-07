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

def generate_random_mutants_separate_files(pssm_df, disorder_df, output_dir="./Fasta", max_mutations=5, mutants_per_size=100, n_attempts=500):
    import os
    from collections import defaultdict

    os.makedirs(output_dir, exist_ok=True)
    mutation_log_path = os.path.join(output_dir, "mutation_log.txt")
    open(mutation_log_path, "w").close()

    if 'SeqNo' in disorder_df.columns:
        disorder_df = disorder_df.rename(columns={'SeqNo': 'Position'})

    merged_df = pd.merge(pssm_df, disorder_df, on='Position')
    disorder_promoting = set(['E', 'K', 'R', 'Q', 'S', 'P', 'G'])
    order_promoting = ['L', 'I', 'F', 'V', 'Y', 'W', 'M', 'C']

    original_sequence = list(merged_df['Residue'])
    sequence_length = len(original_sequence)

    mutation_counts = defaultdict(int)
    mutant_counter = defaultdict(int)

    fasta_files = {
        size: open(os.path.join(output_dir, f"Mutants_{size}_Res.fasta"), "w")
        for size in range(1, max_mutations + 1)
    }

    valid_positions = merged_df[
        (merged_df['Disordered'] == 1) &
        (merged_df['Residue'].isin(disorder_promoting)) &
        (merged_df['DisorderProb'] > 0.7)
    ]['Position'].tolist()

    for block_size in range(1, max_mutations + 1):
        attempts = 0
        while mutation_counts[block_size] < mutants_per_size and attempts < n_attempts:
            attempts += 1
            if len(valid_positions) < block_size:
                break

            selected_positions = random.sample(valid_positions, block_size)
            selected_rows = merged_df[merged_df['Position'].isin(selected_positions)].sort_values('Position')

            mutated_sequence = original_sequence.copy()
            mutation_details = []
            success = True

            for _, row in selected_rows.iterrows():
                pos = row['Position']
                original_aa = row['Residue']
                # scores = row[order_promoting].dropna().astype(float).sort_values(ascending=False)

                # alt_aa = next((aa for aa in scores.index if aa != original_aa), None)


                pssm_threshold = -1.0  # Set your desired threshold

                # Sort alternative residues by score, excluding the original
                scores = row[order_promoting].dropna().astype(float).sort_values(ascending=False)

                # Pick first valid mutation with acceptable PSSM score
                alt_aa = None
                for aa in scores.index:
                    if aa != original_aa and scores[aa] > pssm_threshold:
                        alt_aa = aa
                        break




                if alt_aa:
                    mutated_sequence[pos - 1] = alt_aa
                    mutation_details.append(f"{original_aa}{pos}{alt_aa}")
                else:
                    success = False
                    break

            if not success:
                continue

            mutant_counter[block_size] += 1
            header = f">Mutant_{mutant_counter[block_size]}_Res{block_size}_" + "_".join(mutation_details)

            fasta_file = fasta_files[block_size]
            fasta_file.write(header + "\n")
            for i in range(0, len(mutated_sequence), 60):
                fasta_file.write(''.join(mutated_sequence[i:i + 60]) + "\n")

            with open(mutation_log_path, "a") as log_file:
                log_file.write(f"{header}\n")

            mutation_counts[block_size] += 1

    for f in fasta_files.values():
        f.close()

    print("Random Mutation Summary:")
    for size in range(1, max_mutations + 1):
        print(f" - {mutation_counts[size]} mutants with {size} random mutation(s)")
    print(f"\nFASTA files saved to: {output_dir}")






# === Main Runner ===
if __name__ == "__main__":
    pssm_path = 'fonPB_O.pssm'
    disorder_path = './Dispred/fonPB_O.dispred'
    output_dir = './Fasta'

    pssm_df = read_blast_pssm(pssm_path)
    disorder_df = read_disorder_file(disorder_path)

    generate_random_mutants_separate_files(pssm_df, disorder_df, output_dir=output_dir, max_mutations=5)
