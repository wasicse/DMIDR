# Investigation of intrinsically disordered regions in the Drosophila matrisome

A reusable pipeline for investigating intrinsically disordered regions (IDRs) in the *Drosophila melanogaster* matrisome through sequence analysis, mutant generation, disorder prediction, and optional molecular dynamics (MD) comparison between wild-type and mutant proteins.


## Overview

This pipeline supports the following workflow:

- generate a PSI-BLAST PSSM from a protein sequence
- run DisPredict on the original sequence
- generate consecutive mutants from PSSM and disorder predictions
- run DisPredict on mutant FASTA files
- compare original and mutant disorder predictions
- optionally analyze WT vs mutant MD trajectories

The main goal is to identify candidate sequence regions and genes whose disordered segments may contribute to embryonic matrix dynamics in the *Drosophila* matrisome.

## Project Layout

```text
organized_idr_project/
├── README.md
├── pyproject.toml
├── .python-version
├── configs/
│   └── example.env
├── data/
│   ├── input/
│   │   ├── sequences/
│   │   ├── pdb/
│   │   └── trajectories/
│   └── intermediate/
│       ├── pssm/
│       ├── dispred/
│       └── fasta/
├── results/
├── mdp/
├── notebooks/
├── legacy/
├── scripts/
│   ├── run_pssm.sh
│   ├── run_dispredict.sh
│   ├── run_pipeline.sh
│   ├── generate_mutants.py
│   ├── analyze_disorder.py
│   ├── analyze_md.py
│   └── compare_md_runs.py
└── src/
    └── idr_project/
        ├── io_utils.py
        ├── mutant_generator.py
        ├── disorder_analysis.py
        └── md_analysis.py
```

## Requirements

### System requirements

- Linux or macOS recommended
- Python 3.10
- `psiblast` installed and available in your shell
- access to a BLAST database
- Docker for running DisPredict
- `mkdssp` if you want secondary-structure or ASA calculations

### Python dependencies

This project uses `uv` for Python environment and dependency management.

Install `uv` first, then from the project root run:

```bash
uv sync
```

This will create the environment and install the dependencies listed in `pyproject.toml`.

## Using uv

### Create and sync the environment

```bash
uv sync
```

### Run Python scripts

```bash
uv run python scripts/generate_mutants.py --help
uv run python scripts/analyze_disorder.py --help
uv run python scripts/analyze_md.py --help
```

### Run the pipeline launcher

```bash
uv run bash scripts/run_pipeline.sh configs/local.env
```

### Add a new dependency

```bash
uv add pandas
```

### Add a development dependency

```bash
uv add --dev jupyter ipykernel ruff
```

### Regenerate the lockfile

```bash
uv lock
```

## Input Files

For a new sequence, place your files here.

### Required

- FASTA sequence: `data/input/sequences/SEQ1.fasta`

### Optional for MD analysis

- WT topology: `data/input/pdb/wt.pdb`
- WT trajectory: `data/input/trajectories/wt.xtc`
- mutant topology: `data/input/pdb/mutant.pdb`
- mutant trajectory: `data/input/trajectories/mutant.xtc`

Replace `SEQ1` with your sequence or project name.

## Quick Start

Assume your sequence name is `SEQ1`.

### 1. Generate the PSSM

```bash
bash scripts/run_pssm.sh \
  data/input/sequences/SEQ1.fasta \
  /absolute/path/to/blast/database/nr \
  data/intermediate/pssm/SEQ1.pssm
```

### 2. Run disorder prediction on the original sequence

```bash
docker run -ti -d --name seq1_orig wasicse/dispredict3.0:latest
docker cp data/input/sequences/SEQ1.fasta seq1_orig:/opt/Dispredict3.0/example/sample.fasta
docker exec -i seq1_orig /bin/bash -lc \
  "source /opt/Dispredict3.0/.venv/bin/activate && \
   /opt/Dispredict3.0/.venv/bin/python /opt/Dispredict3.0/script/Dispredict3.0.py \
   -f /opt/Dispredict3.0/example/sample.fasta \
   -o /opt/Dispredict3.0/output/"
docker cp seq1_orig:/opt/Dispredict3.0/output/sample_disPred.txt \
  data/intermediate/dispred/SEQ1_original.dispred
```

### 3. Generate mutant FASTA files

```bash
uv run python scripts/generate_mutants.py \
  --pssm data/intermediate/pssm/SEQ1.pssm \
  --disorder data/intermediate/dispred/SEQ1_original.dispred \
  --output-dir data/intermediate/fasta/SEQ1 \
  --sequence-name SEQ1 \
  --max-mutations 5
```

Typical outputs:

- `data/intermediate/fasta/SEQ1/SEQ1_mutants_1res.fasta`
- `data/intermediate/fasta/SEQ1/SEQ1_mutants_2res.fasta`
- `data/intermediate/fasta/SEQ1/SEQ1_mutants_3res.fasta`
- `data/intermediate/fasta/SEQ1/SEQ1_mutation_log.txt`

### 4. Run disorder prediction on mutant FASTA files

```bash
bash scripts/run_dispredict.sh \
  data/intermediate/fasta/SEQ1 \
  data/intermediate/dispred \
  wasicse/dispredict3.0:latest \
  seq1run \
  5
```

### 5. Compare original vs mutant disorder predictions

```bash
uv run python scripts/analyze_disorder.py \
  --original data/intermediate/dispred/SEQ1_original.dispred \
  --mutant data/intermediate/dispred/SEQ1_mutants_3res.dispred \
  --output-dir results/SEQ1/disorder_3res \
  --label SEQ1_3res
```

### 6. MD analysis

```bash
uv run python scripts/analyze_md.py \
  --wt-top data/input/pdb/wt.pdb \
  --wt-traj data/input/trajectories/wt.xtc \
  --mut-top data/input/pdb/mutant.pdb \
  --mut-traj data/input/trajectories/mutant.xtc \
  --disorder-table results/SEQ1/disorder_probabilities.csv \
  --outdir results/SEQ1/md
```

### 7. Generate a PDB from FASTA (ColabFold)

```bash
bash scripts/run_colabfold.sh \
  data/input/sequences/example.fasta \
  results/example/alphafold
```

You can pass extra ColabFold options after the output directory. The script defaults to:

- input FASTA: /home/mkabir3/Research/47_Dr_Raj_Colab/DMIDR/input/example.fasta
- output dir: result/alphafold

## ASA calculation

Accessible surface area (ASA) is computed per residue using mdtraj's Shrake-Rupley algorithm — no external binary (mkdssp) required. Secondary structure is assigned with mdtraj's built-in DSSP implementation. Relative ASA is $rASA = ASA / MAX\_ASA$ (Wilke scale); residues are labelled "Exposed" if $rASA > 0.25$.

**Single structure:**

```bash
uv run python scripts/calculate_asa.py \
  --pdb data/input/pdb/wt.pdb \
  --out results/SEQ1/wt_asa.csv
```

**All structures in a directory (e.g., ColabFold output):**

```bash
bash scripts/run_asa_all.sh \
  results/example/alphafold \
  results/example/asa
```

Optional chain filter:

```bash
uv run python scripts/calculate_asa.py \
  --pdb data/input/pdb/wt.pdb \
  --out results/SEQ1/wt_chainA_asa.csv \
  --chain A
```

## One-Command Run

1. Copy the example config:

```bash
cp configs/example.env configs/local.env
```

2. Edit `configs/local.env` and set at least:

- `SEQ_NAME`
- `INPUT_FASTA`
- `BLAST_DB`
- `DISPREDICT_IMAGE`
- `MAX_MUTATIONS`

3. Run the workflow:

```bash
uv run bash scripts/run_pipeline.sh configs/local.env
```

### Optional skip flags

You can skip stages by setting any of the following to `1` in the config file:

- `SKIP_PSSM`
- `SKIP_ORIGINAL_DISPREDICT`
- `SKIP_MUTANT_GENERATION`
- `SKIP_MUTANT_DISPREDICT`

## Outputs

Typical outputs include:

### Intermediate

- `.pssm` files
- mutant FASTA files
- original and mutant `.dispred` files

### Final

- disorder probability comparison CSV files
- disorder label comparison CSV files
- summary text files
- disorder comparison plots
- optional MD analysis outputs


## Notes

- Run commands from the project root directory
- Keep all project folders inside the same main directory
- Use sequence-specific names to avoid overwriting previous runs
- If a Docker container with the same name already exists, remove it or use a new name
- The end-to-end workflow depends on your BLAST database path, Docker setup, and available structure or trajectory files


## Citation

If you use this project, please cite:

N Kasirosafar, Md Wasi Ul Kabir, and Md Tamjidul Hoque. *Investigation of intrinsically disordered regions in the Drosophila matrisome*. Wichita State University.


## Authors & Contact

Md Wasi Ul Kabir, Nazanin Kasirosafar, Md Tamjidul Hoque, Raj Logan

Questions/Issues: Md Tamjidul Hoque — thoque@uno.edu
