#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP="$ROOT/tmp/fastcnmf_harmony2"
INPUT_BASE="$TMP/input/gbm_lowres_harmony20_fastcnmf"
COUNTS="$INPUT_BASE.Corrected.HVG.Varnorm.h5ad"
TPM="$INPUT_BASE.TP10K.h5ad"
GENES="$INPUT_BASE.Corrected.HVGs.txt"
OUT="$TMP/parallel"
LOG="$TMP/logs/parallel"
CNMF="$ROOT/.venv-cnmf/bin/cnmf"
TIME="/usr/bin/time"
NAME="gbm_lowres_harmony20_fastcnmf"
WORKERS=4

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

rm -rf "$OUT" "$LOG"
mkdir -p "$OUT" "$LOG"

run_timed() {
  local label="$1"
  shift
  echo "[$(date --iso-8601=seconds)] START $label" | tee "$LOG/${label}.cmd.log"
  printf '%q ' "$@" | tee -a "$LOG/${label}.cmd.log"
  printf '\n' | tee -a "$LOG/${label}.cmd.log"
  "$TIME" -v -o "$LOG/${label}.time.log" "$@" > "$LOG/${label}.stdout.log" 2> "$LOG/${label}.stderr.log"
  echo "[$(date --iso-8601=seconds)] END $label" | tee -a "$LOG/${label}.cmd.log"
}

run_timed prepare "$CNMF" prepare \
  --output-dir "$OUT" --name "$NAME" --counts "$COUNTS" \
  --tpm "$TPM" --genes-file "$GENES" \
  -k 6 8 --n-iter 8 --seed 20260614 --numgenes 3000 --max-nmf-iter 200

echo "[$(date --iso-8601=seconds)] START factorize_${WORKERS}workers" | tee "$LOG/factorize_parallel.cmd.log"
start_ts="$(date +%s)"
status=0
for worker_i in $(seq 0 $((WORKERS - 1))); do
  (
    "$TIME" -v -o "$LOG/factorize_worker_${worker_i}.time.log" \
      "$CNMF" factorize --output-dir "$OUT" --name "$NAME" \
      --total-workers "$WORKERS" --worker-index "$worker_i" \
      > "$LOG/factorize_worker_${worker_i}.stdout.log" \
      2> "$LOG/factorize_worker_${worker_i}.stderr.log"
  ) &
done
for job in $(jobs -p); do
  if ! wait "$job"; then
    status=1
  fi
done
end_ts="$(date +%s)"
{
  echo "workers=$WORKERS"
  echo "elapsed_seconds=$((end_ts - start_ts))"
  echo "exit_status=$status"
  echo "[$(date --iso-8601=seconds)] END factorize_${WORKERS}workers"
} | tee -a "$LOG/factorize_parallel.cmd.log"
if [[ "$status" != "0" ]]; then
  exit "$status"
fi

run_timed combine "$CNMF" combine --output-dir "$OUT" --name "$NAME"
run_timed consensus "$CNMF" consensus --output-dir "$OUT" --name "$NAME"
run_timed k_selection_plot "$CNMF" k_selection_plot --output-dir "$OUT" --name "$NAME"

