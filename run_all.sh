#!/usr/bin/env bash
# Full pipeline runner — processes every sequence in the input FASTA.
#
# Stages (per sequence):
#   1. PSSM              — PSI-BLAST position-specific scoring matrix
#   2. ESMDisPred orig   — disorder prediction on the wild-type sequence
#   3. Mutant generation — consecutive-mutation FASTA files
#   4. ESMDisPred muts   — disorder prediction on each mutant batch
#   5. Disorder analysis — comparison plots and CSV summaries
#   6. ColabFold         — structure prediction (uses pre-computed MSAs if available)
#   7. ASA               — per-residue accessible surface area on ColabFold PDBs
#   8. MD analysis       — RMSF, contacts, DSSP (requires trajectory files)
#
# Usage (from project root):
#   bash run_all.sh [--force] [configs/local.env]
#
# --force   Rerun all selected steps even if outputs already exist.
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"

FORCE=0
CONFIG_FILE="$PROJECT_ROOT/configs/local.env"
for arg in "$@"; do
  if [[ "$arg" == "--force" ]]; then
    FORCE=1
  else
    CONFIG_FILE="$arg"
  fi
done

# ── Interactive step selection ─────────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Pipeline Steps"
echo "================================================================"
echo "  1. PSSM              — PSI-BLAST scoring matrix"
echo "  2. ESMDisPred orig   — disorder prediction (wild-type)"
echo "  3. Mutant generation — consecutive-mutation FASTA files"
echo "  4. ESMDisPred muts   — disorder prediction (mutants)"
echo "  5. Disorder analysis — comparison plots and CSVs"
echo "  6. ColabFold         — structure prediction"
echo "  7. ASA               — accessible surface area"
echo "  8. MD analysis       — RMSF, contacts, DSSP"
echo ""
echo "Enter step numbers to run (e.g. 1 2 3  or  1-5  or  6 7)."
echo "Press Enter with no input to run ALL steps."
echo ""
read -rp "Steps: " STEP_INPUT

# Parse selection into a set of enabled step numbers
declare -A RUN_STEP
for s in 1 2 3 4 5 6 7 8; do RUN_STEP[$s]=0; done

if [[ -z "$STEP_INPUT" ]]; then
  for s in 1 2 3 4 5 6 7 8; do RUN_STEP[$s]=1; done
else
  for token in $STEP_INPUT; do
    if [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      for (( s=${BASH_REMATCH[1]}; s<=${BASH_REMATCH[2]}; s++ )); do
        RUN_STEP[$s]=1
      done
    elif [[ "$token" =~ ^[0-9]+$ ]]; then
      RUN_STEP[$token]=1
    fi
  done
fi

echo ""
echo "Running steps: $(for s in 1 2 3 4 5 6 7 8; do [[ ${RUN_STEP[$s]} -eq 1 ]] && printf '%s ' "$s"; done)"
echo "================================================================"
echo ""

if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
else
  echo "[WARN] Config file not found: $CONFIG_FILE" >&2
  echo "       Copy configs/example.env to configs/local.env and fill in the values." >&2
fi

: "${BLAST_DB:?BLAST_DB is required — set it in $CONFIG_FILE}"
: "${ESMDISPRED_IMAGE:?ESMDISPRED_IMAGE is required — set it in $CONFIG_FILE}"
: "${LARGE_MODELS_DIR:?LARGE_MODELS_DIR is required — set it in $CONFIG_FILE}"

INPUT_FASTA="${INPUT_FASTA:-$PROJECT_ROOT/data/input/sequences/example.fasta}"
MAX_MUTATIONS="${MAX_MUTATIONS:-5}"
MIN_DISORDER_PROB="${MIN_DISORDER_PROB:-0.5}"
REQUIRED_RATIO="${REQUIRED_RATIO:-0.8}"
MODEL="${MODEL:-3}"

# Map menu selection to SKIP_* flags — done after config so menu always wins
SKIP_PSSM=$(( 1 - RUN_STEP[1] ))
SKIP_ORIG_PRED=$(( 1 - RUN_STEP[2] ))
SKIP_MUTANT_GEN=$(( 1 - RUN_STEP[3] ))
SKIP_MUTANT_PRED=$(( 1 - RUN_STEP[4] ))
SKIP_ANALYSIS=$(( 1 - RUN_STEP[5] ))
SKIP_COLABFOLD=$(( 1 - RUN_STEP[6] ))
SKIP_ASA=$(( 1 - RUN_STEP[7] ))
SKIP_MD=$(( 1 - RUN_STEP[8] ))


SEQUENCES_DIR="$PROJECT_ROOT/data/input/sequences"
MSA_DIR="$PROJECT_ROOT/data/input/msas"
PSSM_DIR="$PROJECT_ROOT/data/intermediate/pssm"
DISPRED_DIR="$PROJECT_ROOT/data/intermediate/dispred"

mkdir -p "$SEQUENCES_DIR" "$PSSM_DIR" "$DISPRED_DIR"

# ── ESMDisPred model download (one-time, skipped if already present) ───────────
if [[ ! -f "$LARGE_MODELS_DIR/best.pt" ]]; then
  echo "ESMDisPred models not found — downloading now..."
  bash "$PROJECT_ROOT/scripts/tools/ESMDisPred/run_downloadLargeModels.sh"
fi

# ── Step 0: split multi-FASTA ──────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "Step 0: Splitting $INPUT_FASTA"
echo "================================================================"
uv run --project "$PROJECT_ROOT" python "$PROJECT_ROOT/scripts/split_fasta.py" \
  --input "$INPUT_FASTA" \
  --output-dir "$SEQUENCES_DIR"

# ── Per-sequence loop ──────────────────────────────────────────────────────────
for fasta in "$SEQUENCES_DIR"/*.fasta; do
  [[ "$fasta" -ef "$INPUT_FASTA" ]] && continue   # skip the master input file itself
  seq_name="$(basename "$fasta" .fasta)"
  pssm_out="$PSSM_DIR/${seq_name}.pssm"
  orig_caid="$DISPRED_DIR/${seq_name}_original.caid"
  fasta_dir="$PROJECT_ROOT/data/intermediate/mutants/$seq_name"
  results_dir="$PROJECT_ROOT/results/$seq_name"
  structure_dir="$results_dir/structure"
  asa_dir="$results_dir/asa"
  md_dir="$results_dir/md"

  echo ""
  echo "================================================================"
  echo "Sequence: $seq_name"
  echo "================================================================"
  mkdir -p "$fasta_dir" "$results_dir"

  # ── 1. PSSM ────────────────────────────────────────────────────────────────
  if [[ "$SKIP_PSSM" == "1" ]]; then
    echo "[1/8] Skipping PSSM"
  elif [[ "$FORCE" == "0" && -f "$pssm_out" ]]; then
    echo "[1/8] PSSM already exists, skipping → $pssm_out"
  else
    echo "[1/8] Generating PSSM → $pssm_out"
    bash "$PROJECT_ROOT/scripts/run_pssm.sh" "$fasta" "$BLAST_DB" "$pssm_out"
  fi

  # ── 2. ESMDisPred on original sequence ─────────────────────────────────────
  if [[ "$SKIP_ORIG_PRED" == "1" ]]; then
    echo "[2/8] Skipping ESMDisPred (original)"
  elif [[ "$FORCE" == "0" && -f "$orig_caid" ]]; then
    echo "[2/8] ESMDisPred (original) already exists, skipping → $orig_caid"
  else
    echo "[2/8] ESMDisPred (original) → $orig_caid"
    bash "$PROJECT_ROOT/scripts/run_esmdispred_single.sh" \
      "$fasta" "$DISPRED_DIR" "$ESMDISPRED_IMAGE" "$LARGE_MODELS_DIR" \
      "${seq_name}_original" "$MODEL"
  fi

  # ── 3. Generate mutant FASTA files ─────────────────────────────────────────
  if [[ "$SKIP_MUTANT_GEN" == "1" ]]; then
    echo "[3/8] Skipping mutant generation"
  elif [[ "$FORCE" == "0" ]] && compgen -G "$fasta_dir/*res.fasta" > /dev/null 2>&1; then
    echo "[3/8] Mutant FASTAs already exist, skipping → $fasta_dir"
  else
    echo "[3/8] Generating mutants → $fasta_dir"
    uv run --project "$PROJECT_ROOT" python "$PROJECT_ROOT/scripts/generate_mutants.py" \
      --pssm "$pssm_out" \
      --disorder "$orig_caid" \
      --output-dir "$fasta_dir" \
      --sequence-name "$seq_name" \
      --max-mutations "$MAX_MUTATIONS" \
      --min-disorder-prob "$MIN_DISORDER_PROB" \
      --required-ratio "$REQUIRED_RATIO"
  fi

  # ── 4. ESMDisPred on mutant FASTAs ─────────────────────────────────────────
  if [[ "$SKIP_MUTANT_PRED" == "1" ]]; then
    echo "[4/8] Skipping ESMDisPred (mutants)"
  elif [[ "$FORCE" == "0" && -f "$DISPRED_DIR/${seq_name}_mutants_1res.caid" ]]; then
    echo "[4/8] Mutant predictions already exist, skipping → $DISPRED_DIR"
  else
    echo "[4/8] ESMDisPred (mutants) → $DISPRED_DIR"
    bash "$PROJECT_ROOT/scripts/run_esmdispred.sh" \
      "$fasta_dir" "$DISPRED_DIR" "$ESMDISPRED_IMAGE" "$LARGE_MODELS_DIR" \
      "$seq_name" "$MAX_MUTATIONS" "$MODEL"
  fi

  # ── 5. Disorder comparison analysis ────────────────────────────────────────
  if [[ "$SKIP_ANALYSIS" == "1" ]]; then
    echo "[5/8] Skipping disorder analysis"
  else
    echo "[5/8] Disorder analysis → $results_dir"
    for i in $(seq 1 "$MAX_MUTATIONS"); do
      mutant_caid="$DISPRED_DIR/${seq_name}_mutants_${i}res.caid"
      out_dir="$results_dir/disorder_${i}res"
      if [[ ! -f "$mutant_caid" ]]; then
        echo "  No mutant prediction for block size ${i}, skipping"
        continue
      elif [[ "$FORCE" == "0" ]] && compgen -G "$out_dir/*.csv" > /dev/null 2>&1; then
        echo "  Block size ${i} analysis already exists, skipping"
        continue
      fi
      uv run --project "$PROJECT_ROOT" python "$PROJECT_ROOT/scripts/analyze_disorder.py" \
        --original "$orig_caid" \
        --mutant "$mutant_caid" \
        --output-dir "$out_dir" \
        --label "${seq_name}_${i}res"
    done
  fi

  # ── 6. Structure prediction (ColabFold) ────────────────────────────────────
  if [[ "$SKIP_COLABFOLD" == "1" ]]; then
    echo "[6/8] Skipping ColabFold"
  elif [[ "$FORCE" == "0" ]] && compgen -G "$structure_dir/*.pdb" > /dev/null 2>&1; then
    echo "[6/8] ColabFold output already exists, skipping → $structure_dir"
  else
    echo "[6/8] ColabFold → $structure_dir"
    msa_file="$MSA_DIR/${seq_name}.a3m"
    if [[ -f "$msa_file" ]]; then
      tmp_msa="$MSA_DIR/.tmp_${seq_name}_$$"
      mkdir -p "$tmp_msa"
      cp "$msa_file" "$tmp_msa/"
      bash "$PROJECT_ROOT/scripts/run_colabfold.sh" "$tmp_msa" "$structure_dir"
      rm -rf "$tmp_msa"
    else
      bash "$PROJECT_ROOT/scripts/run_colabfold.sh" "$fasta" "$structure_dir"
    fi
  fi

  # ── 7. ASA calculation on ColabFold PDB outputs ────────────────────────────
  if [[ "$SKIP_ASA" == "1" ]]; then
    echo "[7/8] Skipping ASA"
  else
    echo "[7/8] ASA calculation → $asa_dir"
    if compgen -G "$structure_dir/*.pdb" > /dev/null 2>&1; then
      bash "$PROJECT_ROOT/scripts/run_asa_all.sh" "$structure_dir" "$asa_dir"
    else
      echo "  [WARN] No PDB files in $structure_dir — run ColabFold first (step 6)."
    fi
  fi

  # ── 8. MD analysis ─────────────────────────────────────────────────────────
  if [[ "$SKIP_MD" == "1" ]]; then
    echo "[8/8] Skipping MD analysis"
  else
    echo "[8/8] MD analysis → $md_dir"
    wt_top="$PROJECT_ROOT/data/input/pdb/${seq_name}_wt.pdb"
    wt_traj="$PROJECT_ROOT/data/input/trajectories/${seq_name}_wt.xtc"
    mut_top="$PROJECT_ROOT/data/input/pdb/${seq_name}_mutant.pdb"
    mut_traj="$PROJECT_ROOT/data/input/trajectories/${seq_name}_mutant.xtc"
    disorder_table="$results_dir/disorder_1res/${seq_name}_1res_disorder_probability_comparison.csv"

    missing=()
    [[ -f "$wt_top" ]]        || missing+=("$wt_top")
    [[ -f "$wt_traj" ]]       || missing+=("$wt_traj")
    [[ -f "$mut_top" ]]       || missing+=("$mut_top")
    [[ -f "$mut_traj" ]]      || missing+=("$mut_traj")
    [[ -f "$disorder_table" ]] || missing+=("$disorder_table")

    if [[ ${#missing[@]} -gt 0 ]]; then
      echo "  [WARN] MD analysis skipped for $seq_name — missing files:"
      for f in "${missing[@]}"; do echo "    $f"; done
    else
      uv run --project "$PROJECT_ROOT" python "$PROJECT_ROOT/scripts/analyze_md.py" \
        --wt-top "$wt_top" \
        --wt-traj "$wt_traj" \
        --mut-top "$mut_top" \
        --mut-traj "$mut_traj" \
        --disorder-table "$disorder_table" \
        --outdir "$md_dir"
    fi
  fi

  echo "Finished: $seq_name  →  $results_dir"
done

echo ""
echo "================================================================"
echo "All sequences complete."
echo "Results: $PROJECT_ROOT/results/"
echo "================================================================"
