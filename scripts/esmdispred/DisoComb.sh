#!/bin/bash
# Modified DisoComb.sh — full parallelization of all CPU-bound steps.
# Drop-in replacement for the original; mounted into Docker by run_esmdispred*.sh.
#
# Changes vs original:
#   - DisoRDPbind, fMoRFpred, DFLpred run concurrently (they are independent)
#   - PSI-BLAST runs in parallel (xargs -P nproc) with persistent PSSM cache
#   - All per-sequence steps (PSSM cleanup, IUPred, feature gen, logitReg,
#     NNpackage) run in a single parallel xargs loop

set -uo pipefail

DIRNAME=$1

SCRIPT=$(readlink -f "$0")
SHPATH=$(dirname "$SCRIPT")

if [ $# -eq 0 ]; then
    echo "Usage: ./DisoComb.sh fastafile"
    exit
fi

if [ ! -e "$DIRNAME" ]; then
    echo "File doesn't exist: $DIRNAME"
    exit 1
fi

if [ "${DIRNAME:0:1}" = "/" ]; then
    ABSDIR=$(dirname "$DIRNAME")
else
    ABSDIR="$(pwd)/$(dirname "$DIRNAME")"
fi

FILENAME=${DIRNAME##*/}
FNAME=${FILENAME%.*}
TMPDIR=$ABSDIR/"tmp_"$FNAME
SCOREDIR=$ABSDIR/"features_"$FNAME
PREDIR=$ABSDIR/"pred_"$FNAME

mkdir -p "$TMPDIR" "$SCOREDIR" "$PREDIR"

# ── Step 1: Run 3 whole-fasta predictors concurrently ─────────────────────────
(cd "$SHPATH/programs/DisoRDPbind/" && \
    ./DisoRDPbind "$ABSDIR/$FILENAME" "$TMPDIR/${FNAME}_disordpbind.predictions") &
PID_RDP=$!

(cd "$SHPATH/programs/fMoRFpred/" && \
    ./fMoRFpred.sh "$ABSDIR/$FILENAME" "$TMPDIR/${FNAME}_fmorfpred.predictions") &
PID_FMORF=$!

(cd "$SHPATH/programs/DFLpred/" && \
    java -jar DFLpred.jar "$ABSDIR/$FILENAME" "$TMPDIR/${FNAME}_dflpred.predictions") &
PID_DFL=$!

wait $PID_RDP $PID_FMORF $PID_DFL || true

# Build ID list from .seq files created by the predictors
ls "$TMPDIR"/*.seq | while read line; do
    a=${line##*/}
    echo "${a%%.seq}"
done > "$TMPDIR/idlist"

# ── Step 2: Parallel PSI-BLAST with PSSM cache ────────────────────────────────
NPROC=$(nproc)
CACHE_DIR="${PSSM_CACHE:-}"

run_one_blast() {
    local id="$1"
    local TMPDIR="$2"
    local PREDIR="$3"
    local SHPATH="$4"
    local CACHE_DIR="$5"

    local cache_key cached_pssm=""
    cache_key=$(md5sum "$TMPDIR/$id.seq" | awk '{print $1}')
    [[ -n "$CACHE_DIR" ]] && cached_pssm="$CACHE_DIR/${cache_key}.pssm"

    if [[ -n "$cached_pssm" && -f "$cached_pssm" ]]; then
        cp "$cached_pssm" "$TMPDIR/$id.pssm"
        return 0
    fi

    cd "$SHPATH/programs/blast-2.2.24"
    ./bin/psiblast \
        -query "$TMPDIR/$id.seq" \
        -db ./db/swissprot \
        -num_iterations 3 \
        -out "$TMPDIR/$id.out" \
        -out_ascii_pssm "$TMPDIR/$id.pssm" \
        2>/dev/null || true

    if [[ ! -f "$TMPDIR/$id.pssm" ]]; then
        touch "$PREDIR/use_default_pssm_$id"
        "$SHPATH/programs/create_default_pssm" "$TMPDIR/$id.seq" > "$TMPDIR/$id.pssm"
    elif [[ -n "$cached_pssm" ]]; then
        cp "$TMPDIR/$id.pssm" "${cached_pssm}.tmp.$$"
        mv "${cached_pssm}.tmp.$$" "$cached_pssm"
    fi
}

export -f run_one_blast
export PSSM_CACHE="${PSSM_CACHE:-}" TMPDIR PREDIR SHPATH

xargs -P "$NPROC" -I{} bash -c \
    'run_one_blast "$@"' _ {} "$TMPDIR" "$PREDIR" "$SHPATH" "$CACHE_DIR" \
    < "$TMPDIR/idlist"

# ── Step 3: All per-sequence steps in one parallel loop ───────────────────────
OUTDIR="$SHPATH/output"
mkdir -p "$OUTDIR"
export IUPred_PATH="$SHPATH/programs/iupred"

run_one_sequence() {
    local id="$1"
    local TMPDIR="$2"
    local SCOREDIR="$3"
    local PREDIR="$4"
    local SHPATH="$5"
    local FNAME="$6"
    local OUTDIR="$7"

    # PSSM post-processing
    local pssm="$TMPDIR/$id.pssm"
    if [[ -f "$pssm" ]]; then
        sed '1,3d' "$pssm" | tac | sed '1,6d' | tac > "${pssm}.clean"
        mv "${pssm}.clean" "$pssm"
    fi

    # IUPred
    export IUPred_PATH="$SHPATH/programs/iupred"
    "$SHPATH/programs/iupred/iupred" "$TMPDIR/$id.seq" long  | grep -v '^#' > "$TMPDIR/$id.long"
    "$SHPATH/programs/iupred/iupred" "$TMPDIR/$id.seq" short | grep -v '^#' > "$TMPDIR/$id.short"

    # Feature generation
    grep ">$id" -A 2 "$TMPDIR/${FNAME}_dflpred.predictions"    | sed 's/,/ /g' > "$TMPDIR/$id.dfl"
    grep ">$id" -A 7 "$TMPDIR/${FNAME}_disordpbind.predictions" | sed 's/,/ /g' | sed "s/.*binding.*://g" > "$TMPDIR/$id.rdp"
    grep ">$id" -A 4 "$TMPDIR/${FNAME}_fmorfpred.predictions"   | sed 's/,/ /g' > "$TMPDIR/$id.fmorf"

    # logitReg
    (cd "$SHPATH/programs/logReg/" && ./logitReg "$TMPDIR" "$id") || true
    [[ -f "$TMPDIR/$id.score" ]] || return 0
    mv "$TMPDIR/$id.score"    "$SCOREDIR/$id.score"
    [[ -f "$TMPDIR/$id.log.pred" ]] && mv "$TMPDIR/$id.log.pred" "$PREDIR/$id.log.pred" || true

    # NNpackage — per-sequence tmp dir avoids file conflicts between parallel workers
    local nn_tmp
    nn_tmp=$(mktemp -d)
    cut -d $'\t' -f 3-317 "$SCOREDIR/$id.score" > "$nn_tmp/$id.ttscore"
    cut -d $'\t' -f 1-2   "$SCOREDIR/$id.score" > "$nn_tmp/$id.ttindex"
    (cd "$SHPATH/programs/NNpackage/" && python3 Disnet.py "$nn_tmp/$id.ttscore") > "$nn_tmp/$id.ttpreds"
    paste "$nn_tmp/$id.ttindex" "$nn_tmp/$id.ttpreds" > "$PREDIR/$id.nn.pred"
    cp "$nn_tmp"/* "$OUTDIR/"
    rm -rf "$nn_tmp"
}

export -f run_one_sequence
export SHPATH TMPDIR SCOREDIR PREDIR FNAME OUTDIR

xargs -P "$NPROC" -I{} bash -c \
    'run_one_sequence "$@"' _ {} "$TMPDIR" "$SCOREDIR" "$PREDIR" "$SHPATH" "$FNAME" "$OUTDIR" \
    < "$TMPDIR/idlist"

# ── Cleanup ───────────────────────────────────────────────────────────────────
[[ -d "$TMPDIR"   ]] && rm -rf "$TMPDIR"
[[ -d "$SCOREDIR" ]] && rm -rf "$SCOREDIR"
