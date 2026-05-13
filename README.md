# Investigation of intrinsically disordered regions in the Drosophila matrisome

A pipeline for investigating intrinsically disordered regions (IDRs) in the *Drosophila melanogaster* matrisome. For each input sequence it generates PSSM-guided mutants, predicts disorder with ESMDisPred, predicts structure with ColabFold, computes accessible surface area, and optionally runs MD trajectory analysis.

The goal is to identify candidate sequence regions whose disordered segments may contribute to embryonic matrix dynamics in the *Drosophila* matrisome.

---

## Pipeline stages

Each sequence in `data/input/sequences/example.fasta` is processed through eight stages in order:

| Stage | Description |
|---|---|
| 1 | **PSSM** — PSI-BLAST position-specific scoring matrix |
| 2 | **ESMDisPred (original)** — disorder prediction on the wild-type sequence |
| 3 | **Mutant generation** — consecutive-mutation FASTA files guided by PSSM and disorder scores |
| 4 | **ESMDisPred (mutants)** — disorder prediction on each mutant batch |
| 5 | **Disorder analysis** — comparison plots and CSV summaries (original vs mutant) |
| 6 | **ColabFold** — structure prediction using pre-computed MSAs |
| 7 | **ASA** — per-residue accessible surface area on ColabFold PDB outputs |
| 8 | **MD analysis** — RMSF, CA-CA contacts, and DSSP from trajectory files |

---

## Setup

### Requirements

- Linux or macOS
- Python 3.11
- `psiblast` available in your shell
- Docker (for ESMDisPred)
- BLAST database (e.g. NCBI `nr`)
- ESMDisPred Docker image and large model weights

### Install uv and sync the environment

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env   # or restart your shell
uv sync                        # from the project root
```

---

## Running the pipeline

### 1. Configure

Edit `configs/local.env` and set the three required values:

```bash
BLAST_DB=/path/to/blast/nr
ESMDISPRED_IMAGE=wasicse/esmdispred:latest
LARGE_MODELS_DIR=/path/to/esmdispred/largeModels
```

The input FASTA is pre-set to:
```
/home/mkabir3/Research/47_Dr_Raj_Colab/DMIDR/data/input/sequences/example.fasta
```
It contains five sequences: `vkg-PA`, `fon-PA`, `fon-PB`, `fon-PC`, `fon-PD`.

### 2. Run

```bash
bash run_all.sh configs/local.env
```

The script splits the multi-FASTA and runs all eight stages for each sequence.

### Skip flags

To skip a stage during reruns, set the corresponding flag to `1` in `configs/local.env`:

| Flag | Stage skipped |
|---|---|
| `SKIP_PSSM=1` | PSI-BLAST PSSM |
| `SKIP_ORIG_PRED=1` | ESMDisPred on original |
| `SKIP_MUTANT_GEN=1` | Mutant generation |
| `SKIP_MUTANT_PRED=1` | ESMDisPred on mutants |
| `SKIP_ANALYSIS=1` | Disorder comparison |
| `SKIP_COLABFOLD=1` | ColabFold structure prediction |
| `SKIP_ASA=1` | ASA calculation |
| `SKIP_MD=1` | MD trajectory analysis |

---

## Input files

### Always required

| File | Description |
|---|---|
| `data/input/sequences/example.fasta` | Multi-sequence FASTA (all sequences to process) |
| `data/input/msas/<seq>.a3m` | Pre-computed MSA per sequence for ColabFold |

### Required for MD analysis (stage 8)

Place these files before running. `<seq>` is the sequence name (e.g. `fon-PA`):

| File | Description |
|---|---|
| `data/input/pdb/<seq>_wt.pdb` | Wild-type topology |
| `data/input/pdb/<seq>_mutant.pdb` | Mutant topology |
| `data/input/trajectories/<seq>_wt.xtc` | Wild-type trajectory |
| `data/input/trajectories/<seq>_mutant.xtc` | Mutant trajectory |

If any of these are missing for a sequence, stage 8 prints a warning and continues to the next sequence.

### Generating MSAs on a laptop (if the server is behind a firewall)

The ColabFold MSA server requires internet access. If the server network blocks it, generate MSAs on your laptop first:

```bash
# On your laptop
curl -LsSf https://astral.sh/uv/install.sh | sh
pip install "colabfold[alphafold]"
bash scripts/generate_msa_laptop.sh data/input/sequences/example.fasta msas/

# Upload to the server
rsync -avz msas/ <user>@<server>:<project_path>/data/input/msas/
```

---

## Outputs

After a full run, outputs per sequence are organized as follows:

```
data/intermediate/
├── pssm/<seq>.pssm
├── dispred/
│   ├── <seq>_original.caid
│   └── <seq>_mutants_Nres.caid          (one per block size 1..MAX_MUTATIONS)
└── fasta/<seq>/
    ├── <seq>_mutants_1res.fasta
    ├── <seq>_mutants_2res.fasta
    └── <seq>_mutation_log.txt

results/<seq>/
├── disorder_Nres/
│   ├── <label>_disorder_probability_comparison.csv
│   ├── <label>_disorder_label_comparison.csv
│   ├── <label>_disorder_score_plot.png
│   └── <label>_summary.txt
├── alphafold/                            (ColabFold PDBs and scores)
├── asa/                                  (per-residue ASA CSVs)
└── md/                                   (RMSF, contacts, DSSP plots and CSVs)
```

---

## Project layout

```
DMIDR/
├── run_all.sh                       ← entry point
├── configs/
│   ├── example.env                  ← template
│   └── local.env                    ← fill in and use this
├── data/
│   ├── input/
│   │   ├── sequences/
│   │   │   ├── example.fasta        ← multi-sequence input
│   │   │   └── individual/          ← auto-created by run_all.sh
│   │   ├── msas/                    ← pre-computed .a3m files
│   │   ├── pdb/                     ← wt and mutant topology files
│   │   └── trajectories/            ← wt and mutant .xtc files
│   └── intermediate/
│       ├── pssm/
│       ├── dispred/
│       └── fasta/
├── results/
├── scripts/
│   ├── split_fasta.py
│   ├── run_pssm.sh
│   ├── run_esmdispred_single.sh
│   ├── run_esmdispred.sh
│   ├── run_colabfold.sh
│   ├── run_asa_all.sh
│   ├── generate_msa_laptop.sh
│   ├── generate_mutants.py
│   ├── analyze_disorder.py
│   ├── analyze_md.py
│   ├── calculate_asa.py
│   └── compare_md_runs.py
├── src/
│   └── idr_project/
│       ├── io_utils.py
│       ├── mutant_generator.py
│       ├── disorder_analysis.py
│       └── md_analysis.py
├── notebooks/
├── mdp/
└── legacy/
```

---

## Citation

N Kasirosafar, Md Wasi Ul Kabir, and Md Tamjidul Hoque. *Investigation of intrinsically disordered regions in the Drosophila matrisome*. Wichita State University.

## Authors & Contact

Md Wasi Ul Kabir, Nazanin Kasirosafar, Md Tamjidul Hoque, Raj Logan

Questions/Issues: Md Tamjidul Hoque — thoque@uno.edu
