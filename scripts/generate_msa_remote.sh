#!/usr/bin/env bash
# Generate MSAs via a remote SLURM cluster (e.g. jupyterhub.cs.dmz) that can
# reach the ColabFold MSA server, then sync .a3m files back to this server.
#
# Workflow:
#   1. Copy FASTA to remote via scp
#   2. Submit a SLURM job that installs colabfold (once) and runs --msa-only
#   3. Poll until the job finishes
#   4. Rsync .a3m files back to local_msa_dir
#
# Usage:
#   bash scripts/generate_msa_remote.sh <input.fasta> [local_msa_dir]
#
# Required in configs/local.env:
#   MSA_REMOTE_HOST   — e.g. mkabir3@jupyterhub.cs.dmz
#   MSA_REMOTE_DIR    — scratch dir on remote, e.g. /tmp/colabfold_msa or ~/colabfold_scratch
#
# Optional in configs/local.env:
#   MSA_SLURM_PARTITION  — SLURM partition (default: blank = cluster default)
#   MSA_SLURM_TIME       — wall time limit (default: 02:00:00)
#   MSA_SLURM_MEM        — memory request  (default: 8G)
#   MSA_SLURM_CPUS       — CPUs per task   (default: 4)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_FILE="${CONFIG_FILE:-$PROJECT_ROOT/configs/local.env}"
if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

# ── Args ──────────────────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
  echo "Usage: bash scripts/generate_msa_remote.sh <input.fasta> [local_msa_dir]"
  echo ""
  echo "Required in configs/local.env:"
  echo "  MSA_REMOTE_HOST=user@jupyterhub.cs.dmz"
  echo "  MSA_REMOTE_DIR=/scratch/colabfold_msa   # writable dir on remote"
  exit 1
fi

FASTA="$(realpath "$1")"
LOCAL_MSA_DIR="${2:-$PROJECT_ROOT/data/msas}"
FASTA_BASENAME="$(basename "$FASTA")"

: "${MSA_REMOTE_HOST:?MSA_REMOTE_HOST is required — set it in $CONFIG_FILE}"
: "${MSA_REMOTE_DIR:?MSA_REMOTE_DIR is required — set it in $CONFIG_FILE}"

# SLURM resource defaults (override in local.env)
SLURM_PARTITION="${MSA_SLURM_PARTITION:-}"
SLURM_TIME="${MSA_SLURM_TIME:-02:00:00}"
SLURM_MEM="${MSA_SLURM_MEM:-8G}"
SLURM_CPUS="${MSA_SLURM_CPUS:-4}"

REMOTE_FASTA="$MSA_REMOTE_DIR/$FASTA_BASENAME"
REMOTE_OUT="$MSA_REMOTE_DIR/output"
REMOTE_VENV="$HOME/colabfold_venv"   # on remote
REMOTE_LOG="$MSA_REMOTE_DIR/slurm_msa.log"

echo "================================================================"
echo " Remote MSA generation via SLURM"
echo "  Remote host      : $MSA_REMOTE_HOST"
echo "  Remote scratch   : $MSA_REMOTE_DIR"
echo "  Input FASTA      : $FASTA"
echo "  Local MSA dir    : $LOCAL_MSA_DIR"
echo "  SLURM time/mem   : $SLURM_TIME / $SLURM_MEM"
echo "================================================================"

# ── 1. Check SSH ──────────────────────────────────────────────────────────────
echo ""
echo "[1/5] Checking SSH connection ..."
if ! ssh -o BatchMode=yes -o ConnectTimeout=10 "$MSA_REMOTE_HOST" "echo ok" >/dev/null 2>&1; then
  echo "ERROR: Cannot SSH to $MSA_REMOTE_HOST" >&2
  exit 1
fi
echo "  OK"

# ── 2. Copy FASTA to remote ───────────────────────────────────────────────────
echo ""
echo "[2/5] Copying FASTA to remote ..."
ssh "$MSA_REMOTE_HOST" "mkdir -p '$MSA_REMOTE_DIR' '$REMOTE_OUT'"
scp -q "$FASTA" "${MSA_REMOTE_HOST}:${REMOTE_FASTA}"
echo "  Uploaded: $REMOTE_FASTA"

# ── 3. Build and submit SLURM job ─────────────────────────────────────────────
echo ""
echo "[3/5] Submitting SLURM job ..."

# Build the optional partition line
PARTITION_LINE=""
[[ -n "$SLURM_PARTITION" ]] && PARTITION_LINE="#SBATCH --partition=${SLURM_PARTITION}"

JOB_ID=$(ssh "$MSA_REMOTE_HOST" "
cat > /tmp/colabfold_msa_job.sh << 'SLURM_SCRIPT'
#!/bin/bash
#SBATCH --job-name=colabfold_msa
#SBATCH --output=${REMOTE_LOG}
#SBATCH --time=${SLURM_TIME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
${PARTITION_LINE}

set -euo pipefail
VENV=\"${REMOTE_VENV}\"

echo \"=== Job started: \$(date) ===\"
echo \"=== Node: \$(hostname) ===\"

# Install colabfold into venv if not already done
if [[ ! -x \"\$VENV/bin/colabfold_batch\" ]]; then
  echo \"Installing colabfold into \$VENV ...\"
  python3 -m venv \"\$VENV\"
\"\$VENV/bin/pip\" install --quiet 'colabfold[alphafold]'
  echo \"Install complete.\"
else
  echo \"colabfold already installed at \$VENV/bin/colabfold_batch\"
fi

echo \"Running colabfold_batch --msa-only ...\"
\"\$VENV/bin/colabfold_batch\" --msa-only '${REMOTE_FASTA}' '${REMOTE_OUT}'

echo \"=== Job done: \$(date) ===\"
SLURM_SCRIPT

/usr/local/slurm/bin/sbatch /tmp/colabfold_msa_job.sh 2>&1 | grep -oP '(?<=job )\d+'
")

if [[ -z "$JOB_ID" ]]; then
  echo "ERROR: sbatch did not return a job ID." >&2
  echo "Check the remote for errors." >&2
  exit 1
fi
echo "  Job submitted — ID: $JOB_ID"

# ── 4. Poll until job finishes ────────────────────────────────────────────────
echo ""
echo "[4/5] Waiting for SLURM job $JOB_ID to complete ..."
echo "  (polling every 30 s — Ctrl-C to abort)"

POLL_INTERVAL=30
while true; do
  STATE=$(ssh "$MSA_REMOTE_HOST" "
    /usr/local/slurm/bin/squeue -j '$JOB_ID' -h -o '%T' 2>/dev/null || /usr/local/slurm/bin/sacct -j '$JOB_ID' -n -o State%20 2>/dev/null | head -1 | tr -d ' '
  " 2>/dev/null || true)
  STATE="${STATE:-UNKNOWN}"

  case "$STATE" in
    COMPLETED|COMPLETE)
      echo "  Job $JOB_ID completed."
      break
      ;;
    FAILED|CANCELLED|TIMEOUT|OUT_OF_ME+|NODE_FAIL)
      echo ""
      echo "ERROR: Job $JOB_ID ended with state: $STATE" >&2
      echo "Remote log ($REMOTE_LOG):" >&2
      ssh "$MSA_REMOTE_HOST" "cat '$REMOTE_LOG'" 2>/dev/null >&2 || true
      exit 1
      ;;
    ""|UNKNOWN)
      # squeue returns nothing once job disappears — check sacct
      SACCT_STATE=$(ssh "$MSA_REMOTE_HOST" \
        "/usr/local/slurm/bin/sacct -j '$JOB_ID' -n -o State%20 2>/dev/null | head -1 | tr -d ' '" || true)
      if [[ "$SACCT_STATE" == "COMPLETED" ]]; then
        echo "  Job $JOB_ID completed (via sacct)."
        break
      elif [[ -n "$SACCT_STATE" && "$SACCT_STATE" != "RUNNING" && "$SACCT_STATE" != "PENDING" ]]; then
        echo ""
        echo "ERROR: Job $JOB_ID ended with state: $SACCT_STATE" >&2
        ssh "$MSA_REMOTE_HOST" "cat '$REMOTE_LOG'" 2>/dev/null >&2 || true
        exit 1
      fi
      ;;
  esac

  printf "  State: %-15s  [%s]\r" "$STATE" "$(date +%H:%M:%S)"
  sleep "$POLL_INTERVAL"
done

# Show last few lines of remote log
echo ""
echo "  Remote log (tail):"
ssh "$MSA_REMOTE_HOST" "tail -5 '$REMOTE_LOG'" 2>/dev/null | sed 's/^/    /'

# ── 5. Sync .a3m files back ───────────────────────────────────────────────────
echo ""
echo "[5/5] Syncing .a3m files to $LOCAL_MSA_DIR ..."
mkdir -p "$LOCAL_MSA_DIR"
rsync -avz --include="*.a3m" --exclude="*" \
  "${MSA_REMOTE_HOST}:${REMOTE_OUT}/" "$LOCAL_MSA_DIR/"

A3M_COUNT=$(find "$LOCAL_MSA_DIR" -name "*.a3m" | wc -l)

echo ""
echo "================================================================"
echo " Done!  $A3M_COUNT .a3m file(s) in $LOCAL_MSA_DIR"
echo ""
echo " Run ColabFold with these MSAs:"
echo "   bash scripts/run_colabfold.sh $LOCAL_MSA_DIR <output_dir>"
echo ""
echo " Or via the full pipeline (step 6 auto-detects MSAs):"
echo "   PIPELINE_STEPS=6 bash run_all.sh"
echo "================================================================"
