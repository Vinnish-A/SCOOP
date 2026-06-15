#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$ROOT/tmp/cnmf_spatial_harmony_benchmark"
INPUT_BASE="$TMP_DIR/input/gbm_lowres_visium_3samples_harmony"
COUNTS="$INPUT_BASE.Corrected.HVG.Varnorm.h5ad"
TPM="$INPUT_BASE.TP10K.h5ad"
GENES="$INPUT_BASE.Corrected.HVGs.txt"
CNMF="$ROOT/.venv-cnmf/bin/cnmf"
TIME="/usr/bin/time"
NAME="gbm_lowres_harmony_cnmf"
K_VALUES=(6 8)
N_ITER=8
MAX_NMF_ITER=200
SEED=20260614
PARALLEL_WORKERS=4

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

for required in "$COUNTS" "$TPM" "$GENES"; do
  if [[ ! -f "$required" ]]; then
    echo "Missing required input: $required" >&2
    exit 1
  fi
done

run_timed() {
  local label="$1"
  shift
  local log_dir="$1"
  shift
  mkdir -p "$log_dir"
  echo "[$(date --iso-8601=seconds)] START $label" | tee "$log_dir/${label}.cmd.log"
  printf '%q ' "$@" | tee -a "$log_dir/${label}.cmd.log"
  printf '\n' | tee -a "$log_dir/${label}.cmd.log"
  "$TIME" -v -o "$log_dir/${label}.time.log" "$@" \
    > "$log_dir/${label}.stdout.log" \
    2> "$log_dir/${label}.stderr.log"
  echo "[$(date --iso-8601=seconds)] END $label" | tee -a "$log_dir/${label}.cmd.log"
}

run_one_mode_prepare() {
  local out_dir="$1"
  local log_dir="$2"
  run_timed prepare "$log_dir" "$CNMF" prepare \
    --output-dir "$out_dir" --name "$NAME" --counts "$COUNTS" \
    --tpm "$TPM" --genes-file "$GENES" \
    -k "${K_VALUES[@]}" --n-iter "$N_ITER" --seed "$SEED" \
    --numgenes 3000 --max-nmf-iter "$MAX_NMF_ITER"
}

run_serial() {
  local out_dir="$TMP_DIR/serial"
  local log_dir="$TMP_DIR/logs/serial"
  rm -rf "$out_dir" "$log_dir"
  mkdir -p "$out_dir" "$log_dir"

  run_one_mode_prepare "$out_dir" "$log_dir"
  run_timed factorize "$log_dir" "$CNMF" factorize \
    --output-dir "$out_dir" --name "$NAME" \
    --total-workers 1 --worker-index 0
  run_timed combine "$log_dir" "$CNMF" combine \
    --output-dir "$out_dir" --name "$NAME"
  run_timed consensus "$log_dir" "$CNMF" consensus \
    --output-dir "$out_dir" --name "$NAME"
  run_timed k_selection_plot "$log_dir" "$CNMF" k_selection_plot \
    --output-dir "$out_dir" --name "$NAME"
}

run_parallel() {
  local out_dir="$TMP_DIR/parallel"
  local log_dir="$TMP_DIR/logs/parallel"
  rm -rf "$out_dir" "$log_dir"
  mkdir -p "$out_dir" "$log_dir"

  run_one_mode_prepare "$out_dir" "$log_dir"
  echo "[$(date --iso-8601=seconds)] START factorize_${PARALLEL_WORKERS}workers" \
    | tee "$log_dir/factorize_parallel.cmd.log"
  local start_ts end_ts status
  start_ts="$(date +%s)"
  status=0
  for worker_i in $(seq 0 $((PARALLEL_WORKERS - 1))); do
    (
      "$TIME" -v -o "$log_dir/factorize_worker_${worker_i}.time.log" \
        "$CNMF" factorize \
          --output-dir "$out_dir" --name "$NAME" \
          --total-workers "$PARALLEL_WORKERS" --worker-index "$worker_i" \
        > "$log_dir/factorize_worker_${worker_i}.stdout.log" \
        2> "$log_dir/factorize_worker_${worker_i}.stderr.log"
    ) &
  done
  for job in $(jobs -p); do
    if ! wait "$job"; then
      status=1
    fi
  done
  end_ts="$(date +%s)"
  {
    echo "workers=$PARALLEL_WORKERS"
    echo "elapsed_seconds=$((end_ts - start_ts))"
    echo "exit_status=$status"
    echo "[$(date --iso-8601=seconds)] END factorize_${PARALLEL_WORKERS}workers"
  } | tee -a "$log_dir/factorize_parallel.cmd.log"
  if [[ "$status" != "0" ]]; then
    exit "$status"
  fi

  run_timed combine "$log_dir" "$CNMF" combine \
    --output-dir "$out_dir" --name "$NAME"
  run_timed consensus "$log_dir" "$CNMF" consensus \
    --output-dir "$out_dir" --name "$NAME"
  run_timed k_selection_plot "$log_dir" "$CNMF" k_selection_plot \
    --output-dir "$out_dir" --name "$NAME"
}

case "${1:-all}" in
  serial) run_serial ;;
  parallel) run_parallel ;;
  all)
    run_serial
    run_parallel
    ;;
  *)
    echo "Usage: $0 [serial|parallel|all]" >&2
    exit 2
    ;;
esac
