# Investigation of intrinsically disordered regions in the Drosophila matrisome

A pipeline for investigating intrinsically disordered regions (IDRs) in the *Drosophila melanogaster* matrisome. For each input sequence it generates PSSM-guided mutants, predicts disorder with ESMDisPred, selects the best candidates, predicts structure with ColabFold, computes accessible surface area, and optionally runs MD trajectory analysis.

The goal is to identify candidate sequence variants whose disordered segments are most reduced relative to the wild-type, as potential modulators of embryonic matrix dynamics.

---

## Pipeline stages

| Stage | Script | Description |
|---|---|---|
| 1 | `run_pssm.sh` | **PSSM** — PSI-BLAST position-specific scoring matrix |
| 2 | `run_esmdispred.sh` | **ESMDisPred (wild-type)** — disorder prediction on the original sequence |
| 3 | `generate_mutants.py` | **Mutant generation** — random (default) or consecutive mutations guided by PSSM and disorder scores |
| 4 | `run_esmdispred.sh` | **ESMDisPred (mutants)** — disorder prediction on every mutant batch |
| 5 | `analyze_disorder.py` | **Disorder analysis** — comparison plots and CSV summaries (WT vs mutants) |
| 5.5 | `select_candidates.py` | **Candidate selection** — rank all mutants by avg disorder reduction; write top-N FASTA files |
| 6 | `run_colabfold.sh` | **ColabFold** — structure prediction for WT and all selected candidates |
| 7 | `run_asa_all.sh` | **ASA** — per-residue accessible surface area on ColabFold PDB outputs |
| 8 | `analyze_md.py` | **MD analysis** — RMSF, CA-CA contacts, and DSSP from trajectory files |

---

## Setup

### Requirements

- Linux (Ubuntu 20.04+)
- Python 3.11 (pinned via `.python-version`)
- [uv](https://github.com/astral-sh/uv) package manager
- Docker with NVIDIA GPU support (for ESMDisPred)
- NCBI BLAST+ `psiblast` binary
- BLAST database (e.g. NCBI `nr`)
- SSH access to a SLURM cluster that can reach `api.colabfold.com` (for MSA generation)

### 1. Install uv

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env   # or restart shell
```

### 2. Sync the Python environment

```bash
uv sync   # installs into .venv from pyproject.toml, uses Python 3.11
```

### 3. Install ESMDisPred

```bash
# Clone the tool
git clone https://github.com/wasicse/ESMDisPred scripts/tools/ESMDisPred

# Pull the Docker image
docker pull wasicse/esmdispred:latest

# Download large model weights (~22 GB)
bash scripts/download_esmdispred_models.sh
```

---

## Configuration

Copy the template and fill in required values:

```bash
cp configs/example.env configs/local.env
```

### Required fields

| Variable | Example | Description |
|---|---|---|
| `BLAST_DB` | `/data/blast/nr` | Path prefix to BLAST `nr` database |
| `BLAST_BIN` | `/opt/ncbi-blast-2.11.0+/bin` | Directory containing `psiblast` binary |
| `ESMDISPRED_IMAGE` | `wasicse/esmdispred:latest` | Docker image name |
| `LARGE_MODELS_DIR` | `scripts/tools/ESMDisPred/largeModels` | Path to ESMDisPred model weights |

### Optional overrides (defaults shown)

```bash
INPUT_FASTA=data/input/example.fasta

# Mutation strategy
MUTATION_MODE=random          # random (default) or consecutive
TOP_N_CANDIDATES=3            # how many candidates to select in step 5.5
MIN_DISORDER_PROB=0.5         # residue disorder threshold for eligibility

# Random mode
NUM_VARIANTS=200              # mutant sequences per block size
MUTATIONS_PER_SEQ=5           # max simultaneous mutations (generates 5 FASTA files: 1–5 mut)
MIN_SPACING=5                 # minimum residue spacing between mutations
# MUTATION_SEED=42            # set for reproducibility

# Consecutive mode
MAX_MUTATIONS=5               # block sizes 1–N
REQUIRED_RATIO=0.8            # fraction of block residues that must be disorder-eligible

# ESMDisPred model: 1=ESMDisPred-1  2=ESMDisPred-2  3=ESMDisPred-2PDB  4=ESMDisPred-DNN
MODEL=3

# Skip flags — set to 1 to skip a stage on reruns
SKIP_PSSM=0
SKIP_ORIG_PRED=0
SKIP_MUTANT_GEN=0
SKIP_MUTANT_PRED=0
SKIP_ANALYSIS=0
SKIP_COLABFOLD=0
SKIP_ASA=0
SKIP_MD=0
```

### Remote MSA generation (required if local server has no internet to colabfold)

```bash
MSA_REMOTE_HOST=user@jupyterhub.cs.dmz   # SSH host that can reach api.colabfold.com
MSA_REMOTE_DIR=~/colabfold_scratch        # writable scratch directory on remote
MSA_SLURM_PARTITION=H100                  # SLURM partition on remote
# MSA_SLURM_TIME=02:00:00
# MSA_SLURM_MEM=8G
# MSA_SLURM_CPUS=4
```

When `MSA_REMOTE_HOST` is set and no `.a3m` file exists, step 6 automatically submits a SLURM job to the remote host, waits for it to finish, and rsyncs the resulting `.a3m` back before running ColabFold locally.

---

## Running the pipeline

### Full run (all steps)

```bash
bash run_all.sh configs/local.env
```

### Specific steps

```bash
PIPELINE_STEPS="1 2 3" bash run_all.sh configs/local.env
PIPELINE_STEPS="1-5"   bash run_all.sh configs/local.env
PIPELINE_STEPS="6-8"   bash run_all.sh configs/local.env
```

### Force re-run (skip output-exists checks)

```bash
bash run_all.sh --force configs/local.env
```

### Streamlit dashboard

```bash
uv run streamlit run app.py
```

The dashboard provides:
- **Run & Monitor** tab — launch the pipeline, view live logs, output status per sequence
- **GPU & Containers** tab — GPU utilisation and Docker container progress
- **Results** tab — interactive Plotly disorder profiles, candidate rankings table, download buttons for candidate FASTAs and PDB files

---

## Mutation strategies

### Random (default, `MUTATION_MODE=random`)

Selects `NUM_VARIANTS` random multi-site mutations per block size from the top-100 high-disorder residues, enforcing `MIN_SPACING` between chosen sites. Generates `MUTATIONS_PER_SEQ` FASTA files (one per mutation count 1–N).

### Consecutive (`MUTATION_MODE=consecutive`)

Applies a sliding window over every contiguous block of `N` residues (for N = 1 to `MAX_MUTATIONS`) where at least `REQUIRED_RATIO` of residues exceed `MIN_DISORDER_PROB`. Enumerates all viable positions exhaustively.

Both modes produce the same output format — five FASTA files per sequence with matching header convention.

---

## Candidate selection (step 5.5)

After ESMDisPred runs on all mutants, step 5.5 ranks every mutant by **average disorder reduction**:

```
avg_disorder_reduction = mean(DisorderProb_WT − DisorderProb_mut)   across all residues
```

Secondary metric: `disordered_to_ordered` — count of residues where the binary label flips from 1 (disordered) to 0 (ordered). Candidates with a higher `avg_disorder_reduction` show the greatest global suppression of disorder relative to the wild-type.

The top-`TOP_N_CANDIDATES` mutants are written as individual FASTA files to `outputs/<seq>/candidates/` and a full `candidate_scores.csv` ranking is saved alongside them. ColabFold (step 6) then runs for the wild-type and all selected candidates.

---

## MSA generation (step 6)

ColabFold requires multiple sequence alignments (`.a3m` files). The pipeline resolves them in this order:

1. **Pre-existing** `.a3m` in `data/msas/` — used immediately.
2. **Remote SLURM** — if `MSA_REMOTE_HOST` is set, a job is submitted to a cluster that can reach `api.colabfold.com`. The job installs `colabfold[alphafold]` into a venv on first run, then calls `colabfold_batch --msa-only`. Results are rsynced back automatically.
3. **Local ColabFold** — falls back to running ColabFold directly (requires local MSA server access).

---

## Input files

### Always required

| File | Description |
|---|---|
| `data/input/example.fasta` | Multi-sequence FASTA of all sequences to process |

### Required for MD analysis (stage 8)

`<seq>` is the sequence ID from the FASTA header (e.g. `Muc11A-PB`):

| File | Description |
|---|---|
| `data/pdb/<seq>_wt.pdb` | Wild-type topology |
| `data/pdb/<seq>_mutant.pdb` | Mutant topology |
| `data/trajectories/<seq>_wt.xtc` | Wild-type MD trajectory |
| `data/trajectories/<seq>_mutant.xtc` | Mutant MD trajectory |

If any are missing for a sequence, stage 8 prints a warning and skips that sequence.

---

## Outputs

All pipeline outputs land under `outputs/<seq>/`:

```
outputs/<seq>/
├── pssm/                          PSI-BLAST PSSM file
├── dispred/
│   ├── <seq>_original.caid        WT disorder predictions (pos, AA, prob, label)
│   └── <seq>_mutants_Nres.caid    Mutant predictions (multi-sequence, one block per file)
├── mutants/
│   └── <seq>_mutants_Nres.fasta   Generated mutant sequences
├── disorder_Nres/
│   ├── *_disorder_probability_comparison.csv
│   ├── *_disorder_label_comparison.csv
│   ├── *_disorder_score_plot.png
│   └── *_summary.txt
├── candidates/
│   ├── candidate_scores.csv        Full mutant ranking (avg disorder reduction)
│   └── candidate_N_<id>.fasta      Top-N individual candidate FASTAs
├── structure/                     ColabFold PDB and score files (WT + candidates)
├── asa/                           Per-residue ASA CSVs
└── md/                            RMSF, contacts, DSSP outputs
```

---

## Project layout

```
DMIDR/
├── run_all.sh                          pipeline entry point
├── app.py                              Streamlit dashboard
├── configs/
│   ├── example.env                     template (copy to local.env)
│   └── local.env                       machine-specific config (gitignored)
├── data/
│   ├── input/                          input FASTA files
│   ├── msas/                           pre-computed .a3m files (gitignored)
│   ├── pdb/                            topology PDBs for MD (gitignored)
│   └── trajectories/                   .xtc trajectory files for MD (gitignored)
├── outputs/                            all pipeline outputs (gitignored)
├── scripts/
│   ├── split_fasta.py
│   ├── run_pssm.sh
│   ├── run_esmdispred_single.sh
│   ├── run_esmdispred.sh
│   ├── download_esmdispred_models.sh
│   ├── run_colabfold.sh
│   ├── run_asa_all.sh
│   ├── generate_msa_remote.sh          SLURM-based remote MSA generation
│   ├── generate_mutants.py             mutant generation CLI
│   ├── select_candidates.py            candidate ranking and selection
│   ├── analyze_disorder.py
│   ├── analyze_md.py
│   └── calculate_asa.py
├── src/
│   └── idr_project/
│       ├── mutant_generator.py         random + consecutive mutation logic
│       ├── disorder_analysis.py
│       ├── md_analysis.py
│       └── io_utils.py
├── notebooks/
└── legacy/
```

---

## Citation

N Kasirosafar, Md Wasi Ul Kabir, and Md Tamjidul Hoque. *Investigation of intrinsically disordered regions in the Drosophila matrisome*. Wichita State University.

## Authors & Contact

Md Wasi Ul Kabir, Nazanin Kasirosafar, Md Tamjidul Hoque, Raj Logan

Questions / Issues: Md Tamjidul Hoque — thoque@uno.edu
