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
#   bash run_all.sh [--force] [--clean] [configs/local.env]
#
# --force   Rerun all selected steps even if outputs already exist.
# --clean   Delete all previous outputs before running (fresh start).
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"

# ── Logging — tee everything to a timestamped file in outputs/ ────────────────
mkdir -p "$PROJECT_ROOT/outputs"
LOG_FILE="$PROJECT_ROOT/outputs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "Logging to: $LOG_FILE"

FORCE=0
CLEAN=0
CONFIG_FILE="$PROJECT_ROOT/configs/local.env"
for arg in "$@"; do
  if [[ "$arg" == "--force" ]]; then
    FORCE=1
  elif [[ "$arg" == "--clean" ]]; then
    CLEAN=1
  else
    CONFIG_FILE="$arg"
  fi
done

# Source config early so ESMDISPRED_IMAGE is available for container stop.
if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

# ── Clean previous run outputs ────────────────────────────────────────────────
if [[ "$CLEAN" == "1" ]]; then
  echo ""
  echo "================================================================"
  echo " --clean: removing all generated outputs"
  echo "================================================================"
  # Stop any running ESMDisPred containers (by image name and by output-dir mounts)
  if [[ -n "${ESMDISPRED_IMAGE:-}" ]]; then
    running=$(docker ps --filter "ancestor=${ESMDISPRED_IMAGE}" --format "{{.ID}}" 2>/dev/null || true)
    if [[ -n "$running" ]]; then
      echo "  Stopping ESMDisPred containers: $running"
      docker stop $running 2>/dev/null || true
      sleep 2
    fi
  fi
  # Also stop any container whose mounts include our outputs directory
  extra=$(docker ps --format "{{.ID}} {{.Mounts}}" 2>/dev/null \
    | grep "$PROJECT_ROOT" | awk '{print $1}' || true)
  if [[ -n "$extra" ]]; then
    echo "  Stopping containers with project mounts: $extra"
    docker stop $extra 2>/dev/null || true
    sleep 2
  fi
  # Delete sequence outputs (keep pssm_cache for reuse)
  for seq_dir in "$PROJECT_ROOT/outputs"/*/; do
    [[ "$(basename "$seq_dir")" == "pssm_cache" ]] && continue
    [[ "$(basename "$seq_dir")" == "sequences" ]] && continue
    echo "  Removing $seq_dir"
    rm -rf "$seq_dir" 2>/dev/null || true
  done
  # Delete wrong ESM2/Dispredict3.0 mutant feature caches
  find "$PROJECT_ROOT/scripts/tools/ESMDisPred/features/ESM2" \
       -name "Mutant_*.csv" -delete 2>/dev/null || true
  find "$PROJECT_ROOT/scripts/tools/ESMDisPred/features/Dispredict3.0" \
       -name "Mutant_*" -delete 2>/dev/null || true
  # Delete old log files (except the current one)
  find "$PROJECT_ROOT/outputs" -maxdepth 1 -name "run_*.log" \
       ! -name "$(basename "$LOG_FILE")" -delete 2>/dev/null || true
  echo "  Clean complete."
  echo ""
fi

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
if [[ -n "${PIPELINE_STEPS:-}" ]]; then
  STEP_INPUT="$PIPELINE_STEPS"
  echo "Steps (from env): $STEP_INPUT"
elif [[ -t 0 ]]; then
  read -rp "Steps: " STEP_INPUT
else
  echo "ERROR: PIPELINE_STEPS env var must be set when running non-interactively." >&2
  exit 1
fi

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

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "[WARN] Config file not found: $CONFIG_FILE" >&2
  echo "       Copy configs/example.env to configs/local.env and fill in the values." >&2
fi

: "${BLAST_DB:?BLAST_DB is required — set it in $CONFIG_FILE}"
: "${ESMDISPRED_IMAGE:?ESMDISPRED_IMAGE is required — set it in $CONFIG_FILE}"
: "${LARGE_MODELS_DIR:?LARGE_MODELS_DIR is required — set it in $CONFIG_FILE}"
export BLAST_BIN  # make available to background psiblast subprocesses

INPUT_FASTA="${INPUT_FASTA:-$PROJECT_ROOT/data/input/example.fasta}"
MODEL="${MODEL:-3}"

# Mutation mode: random (default) or consecutive
MUTATION_MODE="${MUTATION_MODE:-random}"
MIN_DISORDER_PROB="${MIN_DISORDER_PROB:-0.5}"

# Random-mode parameters
NUM_VARIANTS="${NUM_VARIANTS:-200}"
MUTATIONS_PER_SEQ="${MUTATIONS_PER_SEQ:-5}"
MIN_SPACING="${MIN_SPACING:-5}"
MUTATION_SEED="${MUTATION_SEED:-}"

# Consecutive-mode parameters
MAX_MUTATIONS="${MAX_MUTATIONS:-5}"
REQUIRED_RATIO="${REQUIRED_RATIO:-0.8}"

# Map menu selection to SKIP_* flags — done after config so menu always wins
SKIP_PSSM=$(( 1 - RUN_STEP[1] ))
SKIP_ORIG_PRED=$(( 1 - RUN_STEP[2] ))
SKIP_MUTANT_GEN=$(( 1 - RUN_STEP[3] ))
SKIP_MUTANT_PRED=$(( 1 - RUN_STEP[4] ))
SKIP_ANALYSIS=$(( 1 - RUN_STEP[5] ))
SKIP_COLABFOLD=$(( 1 - RUN_STEP[6] ))
SKIP_ASA=$(( 1 - RUN_STEP[7] ))
SKIP_MD=$(( 1 - RUN_STEP[8] ))


SEQUENCES_DIR="$PROJECT_ROOT/outputs/sequences"
MSA_DIR="$PROJECT_ROOT/data/msas"
PSSM_CACHE_DIR="$PROJECT_ROOT/outputs/pssm_cache"
export PSSM_CACHE_DIR

mkdir -p "$SEQUENCES_DIR" "$PSSM_CACHE_DIR"

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

# ── Step 1: PSSM — run all sequences in parallel ───────────────────────────────
mapfile -t _all_fastas < <(
  for f in "$SEQUENCES_DIR"/*.fasta; do
    [[ "$f" -ef "$INPUT_FASTA" ]] && continue
    echo "$f"
  done
)
_seq_count="${#_all_fastas[@]}"

if [[ "$SKIP_PSSM" == "1" ]]; then
  echo "[1/8] Skipping PSSM (all sequences)"
elif [[ "$_seq_count" -gt 0 ]]; then
  # Divide cores evenly; cap at 16 per job (PSI-BLAST I/O-bound beyond that)
  _threads_each=$(( $(nproc) / _seq_count ))
  _threads_each=$(( _threads_each < 1 ? 1 : _threads_each > 16 ? 16 : _threads_each ))
  echo ""
  echo "================================================================"
  echo "Step 1: PSSM — ${_seq_count} sequences × ${_threads_each} threads in parallel"
  echo "================================================================"
  _pssm_pids=()
  for _fasta in "${_all_fastas[@]}"; do
    _sname="$(basename "$_fasta" .fasta)"
    _pout="$PROJECT_ROOT/outputs/$_sname/pssm/${_sname}.pssm"
    mkdir -p "$(dirname "$_pout")"
    if [[ "$FORCE" == "0" && -f "$_pout" ]]; then
      echo "  [skip] $_sname — already exists"
    else
      echo "  [run]  $_sname → $_pout"
      bash "$PROJECT_ROOT/scripts/run_pssm.sh" \
        "$_fasta" "$BLAST_DB" "$_pout" "$_threads_each" &
      _pssm_pids+=($!)
    fi
  done
  if [[ ${#_pssm_pids[@]} -gt 0 ]]; then
    echo "  Waiting for ${#_pssm_pids[@]} PSSM job(s)..."
    for _pid in "${_pssm_pids[@]}"; do wait "$_pid" || { echo "PSSM job $_pid failed" >&2; exit 1; }; done
    echo "  All PSSM jobs done."
  fi
fi

# ── Per-sequence loop ──────────────────────────────────────────────────────────
for fasta in "${_all_fastas[@]}"; do
  seq_name="$(basename "$fasta" .fasta)"
  results_dir="$PROJECT_ROOT/outputs/$seq_name"
  pssm_dir="$results_dir/pssm"
  dispred_dir="$results_dir/dispred"
  mutants_dir="$results_dir/mutants"
  structure_dir="$results_dir/structure"
  asa_dir="$results_dir/asa"
  md_dir="$results_dir/md"
  pssm_out="$pssm_dir/${seq_name}.pssm"
  orig_caid="$dispred_dir/${seq_name}_original.caid"

  echo ""
  echo "================================================================"
  echo "Sequence: $seq_name"
  echo "================================================================"
  mkdir -p "$pssm_dir" "$dispred_dir" "$mutants_dir" "$structure_dir" "$asa_dir"

  # ── 1. PSSM (already done above in parallel) ───────────────────────────────
  if [[ "$SKIP_PSSM" == "1" ]]; then
    echo "[1/8] Skipping PSSM"
  elif [[ -f "$pssm_out" ]]; then
    echo "[1/8] PSSM done → $pssm_out"
  else
    echo "[1/8] PSSM missing for $seq_name — something went wrong" >&2
    exit 1
  fi

  # ── 2. ESMDisPred on original sequence ─────────────────────────────────────
  if [[ "$SKIP_ORIG_PRED" == "1" ]]; then
    echo "[2/8] Skipping ESMDisPred (original)"
  elif [[ "$FORCE" == "0" && -f "$orig_caid" ]]; then
    echo "[2/8] ESMDisPred (original) already exists, skipping → $orig_caid"
  else
    echo "[2/8] ESMDisPred (original) → $orig_caid"
    bash "$PROJECT_ROOT/scripts/run_esmdispred_single.sh" \
      "$fasta" "$dispred_dir" "$ESMDISPRED_IMAGE" "$LARGE_MODELS_DIR" \
      "${seq_name}_original" "$MODEL"
  fi

  # ── 3. Generate mutant FASTA files ─────────────────────────────────────────
  if [[ "$SKIP_MUTANT_GEN" == "1" ]]; then
    echo "[3/8] Skipping mutant generation"
  elif [[ "$FORCE" == "0" ]] && compgen -G "$mutants_dir/*res.fasta" > /dev/null 2>&1; then
    echo "[3/8] Mutant FASTAs already exist, skipping → $mutants_dir"
  else
    echo "[3/8] Generating mutants (mode=$MUTATION_MODE) → $mutants_dir"
    _seed_arg=()
    [[ -n "$MUTATION_SEED" ]] && _seed_arg=(--seed "$MUTATION_SEED")
    uv run --project "$PROJECT_ROOT" python "$PROJECT_ROOT/scripts/generate_mutants.py" \
      --pssm "$pssm_out" \
      --disorder "$orig_caid" \
      --output-dir "$mutants_dir" \
      --sequence-name "$seq_name" \
      --mode "$MUTATION_MODE" \
      --min-disorder-prob "$MIN_DISORDER_PROB" \
      --num-variants "$NUM_VARIANTS" \
      --mutations-per-seq "$MUTATIONS_PER_SEQ" \
      --min-spacing "$MIN_SPACING" \
      --max-mutations "$MAX_MUTATIONS" \
      --required-ratio "$REQUIRED_RATIO" \
      "${_seed_arg[@]}"
  fi

  # ── 4. ESMDisPred on mutant FASTAs ─────────────────────────────────────────
  if [[ "$SKIP_MUTANT_PRED" == "1" ]]; then
    echo "[4/8] Skipping ESMDisPred (mutants)"
  else
    # Check each block size individually — only skip sizes that already have a .caid
    _need_pred=0
    for _i in $(seq 1 "$MAX_MUTATIONS"); do
      _fasta=$(find "$mutants_dir" -maxdepth 1 -name "*_mutants_${_i}res.fasta" | head -n 1)
      [[ -z "$_fasta" ]] && continue
      if [[ "$FORCE" == "0" && -f "$dispred_dir/${seq_name}_mutants_${_i}res.caid" ]]; then
        echo "[4/8] Block ${_i}res already exists, skipping"
      else
        _need_pred=1
      fi
    done
    if [[ "$_need_pred" == "1" ]]; then
      echo "[4/8] ESMDisPred (mutants) → $dispred_dir"
      bash "$PROJECT_ROOT/scripts/run_esmdispred.sh" \
        "$mutants_dir" "$dispred_dir" "$ESMDISPRED_IMAGE" "$LARGE_MODELS_DIR" \
        "$seq_name" "$MAX_MUTATIONS" "$MODEL"
    else
      echo "[4/8] All mutant predictions already exist → $dispred_dir"
    fi
  fi

  # ── 5. Disorder comparison analysis ────────────────────────────────────────
  if [[ "$SKIP_ANALYSIS" == "1" ]]; then
    echo "[5/8] Skipping disorder analysis"
  else
    echo "[5/8] Disorder analysis → $results_dir"
    for i in $(seq 1 "$MAX_MUTATIONS"); do
      mutant_caid="$dispred_dir/${seq_name}_mutants_${i}res.caid"
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
    wt_top="$PROJECT_ROOT/data/pdb/${seq_name}_wt.pdb"
    wt_traj="$PROJECT_ROOT/data/trajectories/${seq_name}_wt.xtc"
    mut_top="$PROJECT_ROOT/data/pdb/${seq_name}_mutant.pdb"
    mut_traj="$PROJECT_ROOT/data/trajectories/${seq_name}_mutant.xtc"
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
echo "Results: $PROJECT_ROOT/outputs/"
echo "================================================================"
