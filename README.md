# Investigation of intrinsically disordered regions in the Drosophila matrisome

A pipeline for investigating intrinsically disordered regions (IDRs) in the *Drosophila melanogaster* matrisome. For each input sequence it generates PSSM-guided mutants, predicts disorder with ESMDisPred, predicts structure with ColabFold, computes accessible surface area, and optionally runs MD trajectory analysis.

The goal is to identify candidate sequence regions whose disordered segments may contribute to embryonic matrix dynamics in the *Drosophila* matrisome.

---

## Pipeline stages

Each sequence in `data/input/example.fasta` is processed through eight stages:

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
- Python 3.11+
- `psiblast` available in your shell (via conda or system install)
- Docker with NVIDIA GPU support (for ESMDisPred)
- BLAST database (e.g. NCBI `nr`)

### 1. Install uv

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env   # or restart your shell
```

### 2. Sync the Python environment

```bash
uv sync   # from the project root
```

### 3. Install ESMDisPred tool and models

```bash
# Clone the tool into scripts/tools/
git clone https://github.com/wasicse/ESMDisPred scripts/tools/ESMDisPred

# Pull the Docker image
docker pull wasicse/esmdispred:latest

# Download large model weights (~22 GB)
bash scripts/download_esmdispred_models.sh
```

---

## Running the pipeline

### 1. Configure

Copy the template and fill in required values:

```bash
cp configs/example.env configs/local.env
```

Edit `configs/local.env`:

```bash
BLAST_DB=/path/to/blast/nr              # BLAST database prefix
ESMDISPRED_IMAGE=wasicse/esmdispred:latest
LARGE_MODELS_DIR=/path/to/DMIDR/scripts/tools/ESMDisPred/largeModels
INPUT_FASTA=/path/to/DMIDR/data/input/sequences/example.fasta
```

### 2. Run

```bash
bash run_all.sh [configs/local.env]
```

The script shows an interactive step menu — choose which stages to run:

```
================================================================
 Pipeline Steps
================================================================
  1. PSSM              — PSI-BLAST scoring matrix
  2. ESMDisPred orig   — disorder prediction (wild-type)
  3. Mutant generation — consecutive-mutation FASTA files
  4. ESMDisPred muts   — disorder prediction (mutants)
  5. Disorder analysis — comparison plots and CSVs
  6. ColabFold         — structure prediction
  7. ASA               — accessible surface area
  8. MD analysis       — RMSF, contacts, DSSP

Enter step numbers to run (e.g. 1 2 3  or  1-5  or  6 7).
Press Enter with no input to run ALL steps.
```

Use `--force` to rerun stages even if outputs already exist:

```bash
bash run_all.sh --force configs/local.env
```

### Mutant generation tuning

Three parameters in `configs/local.env` control how many mutants are generated:

| Parameter | Default | Effect |
|---|---|---|
| `MAX_MUTATIONS` | `5` | Number of consecutive-residue block sizes (1–N) |
| `MIN_DISORDER_PROB` | `0.5` | Minimum disorder score for a site to be eligible |
| `REQUIRED_RATIO` | `0.8` | Fraction of residues in a block that must be eligible |

Raise `MIN_DISORDER_PROB` and `REQUIRED_RATIO` to generate fewer, higher-confidence mutants.

---

## Input files

### Always required

| File | Description |
|---|---|
| `data/input/example.fasta` | Multi-sequence FASTA — all sequences to process |

### Optional (for ColabFold without MSA server)

The ColabFold MSA server requires internet access. If your server network blocks it, generate MSAs on a laptop first:

```bash
# On your laptop
curl -LsSf https://astral.sh/uv/install.sh | sh && uv sync
bash scripts/generate_msa_laptop.sh data/input/example.fasta msas/

# Upload to the server
rsync -avz msas/ <user>@<server>:<project_path>/data/msas/
```

Place the resulting `.a3m` files in `data/msas/` — the pipeline will use them automatically.

### Required for MD analysis (stage 8)

Place these files before running stage 8. `<seq>` is the sequence name (e.g. `fon-PA`):

| File | Description |
|---|---|
| `data/pdb/<seq>_wt.pdb` | Wild-type topology |
| `data/pdb/<seq>_mutant.pdb` | Mutant topology |
| `data/trajectories/<seq>_wt.xtc` | Wild-type trajectory |
| `data/trajectories/<seq>_mutant.xtc` | Mutant trajectory |

If any are missing for a sequence, stage 8 prints a warning and skips that sequence.

---

## Project layout

```
DMIDR/
├── run_all.sh                       ← pipeline entry point
├── configs/
│   ├── example.env                  ← template (copy to local.env)
│   └── local.env                    ← your machine-specific config (gitignored)
├── data/
│   ├── input/
│   │   └── example.fasta            ← only file you provide
│   ├── msas/                        ← pre-computed .a3m files (gitignored)
│   ├── pdb/                         ← wt/mutant topology files for MD (gitignored)
│   └── trajectories/                ← wt/mutant .xtc files for MD (gitignored)
├── outputs/                         ← all pipeline outputs (gitignored)
│   ├── sequences/                   ← auto-split per-sequence FASTAs
│   └── <seq>/
│       ├── pssm/                    ← PSI-BLAST PSSM
│       ├── dispred/                 ← ESMDisPred .caid files
│       ├── mutants/                 ← generated mutant FASTA files
│       ├── structure/               ← ColabFold PDB and score files
│       ├── asa/                     ← per-residue ASA CSVs
│       ├── disorder_Nres/           ← disorder comparison plots and CSVs
│       └── md/                      ← RMSF, contacts, DSSP outputs
├── scripts/
│   ├── split_fasta.py
│   ├── run_pssm.sh
│   ├── run_esmdispred_single.sh
│   ├── run_esmdispred.sh
│   ├── download_esmdispred_models.sh
│   ├── run_colabfold.sh
│   ├── run_asa_all.sh
│   ├── generate_msa_laptop.sh
│   ├── generate_mutants.py
│   ├── analyze_disorder.py
│   ├── analyze_md.py
│   └── calculate_asa.py
├── src/
│   └── idr_project/
│       ├── io_utils.py
│       ├── mutant_generator.py
│       ├── disorder_analysis.py
│       └── md_analysis.py
├── notebooks/
└── legacy/
```

---

## Citation

N Kasirosafar, Md Wasi Ul Kabir, and Md Tamjidul Hoque. *Investigation of intrinsically disordered regions in the Drosophila matrisome*. Wichita State University.

## Authors & Contact

Md Wasi Ul Kabir, Nazanin Kasirosafar, Md Tamjidul Hoque, Raj Logan

Questions/Issues: Md Tamjidul Hoque — thoque@uno.edu
