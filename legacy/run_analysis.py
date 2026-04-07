import pandas as pd
import matplotlib.pyplot as plt

# === Function to load disorder prediction files ===
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


# === Load files ===
original_disorder_df = read_disorder_file("./Dispred/fonPB_O.dispred")
mutated_disorder_df = read_disorder_file("./Dispred/fonPB_M.dispred")

# === Merge for analysis ===
comparison_df = pd.merge(original_disorder_df, mutated_disorder_df,
                         on='SeqNo', suffixes=('_orig', '_mut'))

# === Part 1: Compare Disorder Probabilities ===
comparison_df['DeltaDisorder'] = comparison_df['DisorderProb_mut'] - comparison_df['DisorderProb_orig']

# Summary statistics for disorder probabilities
num_residues = len(comparison_df)
avg_delta = comparison_df['DeltaDisorder'].mean()
num_increased = (comparison_df['DeltaDisorder'] > 0).sum()
num_decreased = (comparison_df['DeltaDisorder'] < 0).sum()

# Output: Disorder probability comparison
# print("\nDisorder Probability Comparison Per Residue:")
# print(comparison_df[['SeqNo', 'AA_orig', 'AA_mut', 'DisorderProb_orig', 'DisorderProb_mut', 'DeltaDisorder']])

print("\nSummary Statistics (Disorder Probability):")
print(f"Total residues compared: {num_residues}")
print(f"Average disorder change: {avg_delta:.3f}")
print(f"Residues with increased disorder: {num_increased}")
print(f"Residues with decreased disorder: {num_decreased}")

# Optional: save to CSV
comparison_df[['SeqNo', 'AA_orig', 'AA_mut', 'DisorderProb_orig', 'DisorderProb_mut', 'DeltaDisorder']]\
    .to_csv("disorder_probability_comparison.csv", index=False)

# Plotting disorder probabilities
plt.figure(figsize=(10, 4))
plt.plot(comparison_df['SeqNo'], comparison_df['DisorderProb_orig'], label='Original', marker='o')
plt.plot(comparison_df['SeqNo'], comparison_df['DisorderProb_mut'], label='Mutated', marker='x')
plt.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
plt.title("Disorder Score Comparison (Original vs Mutated)")
plt.xlabel("Residue Position")
plt.ylabel("Disorder Probability")
plt.legend()
plt.tight_layout()
plt.show()

# === Part 2: Compare Disorder Labels ===
comparison_df['LabelChange'] = comparison_df['Disordered_mut'] - comparison_df['Disordered_orig']

# Label change counts
total = len(comparison_df)
unchanged = (comparison_df['LabelChange'] == 0).sum()
to_disordered = (comparison_df['LabelChange'] == 1).sum()
to_ordered = (comparison_df['LabelChange'] == -1).sum()

# # Output: Label change summary
# print("\nLabel Change Summary (0 = ordered, 1 = disordered):")
# print(comparison_df[['SeqNo', 'AA_orig', 'AA_mut', 'Disordered_orig', 'Disordered_mut', 'LabelChange']])

print("\nSummary Statistics (Disorder Labels):")
print(f"Total residues:             {total}")
print(f"Unchanged predictions:      {unchanged}")
print(f"Ordered → Disordered:       {to_disordered}")
print(f"Disordered → Ordered:       {to_ordered}")

# Optional: save to CSV
comparison_df[['SeqNo', 'AA_orig', 'AA_mut', 'Disordered_orig', 'Disordered_mut', 'LabelChange']]\
    .to_csv("disorder_label_comparison.csv", index=False)
