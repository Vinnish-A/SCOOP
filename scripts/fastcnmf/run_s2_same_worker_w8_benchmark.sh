#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_s2_same_worker_w8_benchmark.sh --dataset-id ID --input-h5ad PATH [options]

Options:
  --out-root PATH              Default: tmp/fastcnmf_large_benchmark
  --sample-key KEY             Default: sample_id
  --reference-python PATH      Default: ./.venv-cnmf/bin/python
  --candidate-python PATH      Default: ./.venv-cnmf-h2/bin/python
  --workers N                  Default: 8
  --consensus-workers N        Default: 4
  --k-values "6 8 10 12"       Default: "6 8 10 12"
  --n-iter N                   Default: 20
  --ref-max-nmf-iter N         Default: 200
  --fast-max-nmf-iter N        Default: 50
  --n-top-genes N              Default: 3000
  --seed N                     Default: 20260614
  --force                      Re-run stages even when a successful time log exists
USAGE
}

DATASET_ID=""
INPUT_H5AD=""
OUT_ROOT="tmp/fastcnmf_large_benchmark"
SAMPLE_KEY="sample_id"
REFERENCE_PYTHON="./.venv-cnmf/bin/python"
CANDIDATE_PYTHON="./.venv-cnmf-h2/bin/python"
WORKERS=8
CONSENSUS_WORKERS=4
K_VALUES="6 8 10 12"
N_ITER=20
REF_MAX_NMF_ITER=200
FAST_MAX_NMF_ITER=50
N_TOP_GENES=3000
SEED=20260614
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset-id) DATASET_ID="$2"; shift 2 ;;
    --input-h5ad) INPUT_H5AD="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --sample-key) SAMPLE_KEY="$2"; shift 2 ;;
    --reference-python) REFERENCE_PYTHON="$2"; shift 2 ;;
    --candidate-python) CANDIDATE_PYTHON="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --consensus-workers) CONSENSUS_WORKERS="$2"; shift 2 ;;
    --k-values) K_VALUES="$2"; shift 2 ;;
    --n-iter) N_ITER="$2"; shift 2 ;;
    --ref-max-nmf-iter) REF_MAX_NMF_ITER="$2"; shift 2 ;;
    --fast-max-nmf-iter) FAST_MAX_NMF_ITER="$2"; shift 2 ;;
    --n-top-genes) N_TOP_GENES="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$DATASET_ID" || -z "$INPUT_H5AD" ]]; then
  usage >&2
  exit 2
fi

read -r -a K_ARRAY <<< "$K_VALUES"
if [[ ${#K_ARRAY[@]} -eq 0 ]]; then
  echo "--k-values must contain at least one k" >&2
  exit 2
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=""
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"

CNMF_BIN="$(dirname "$REFERENCE_PYTHON")/cnmf"
DATA_ROOT="$OUT_ROOT/$DATASET_ID"
if [[ "$N_ITER" == "8" ]]; then
  ITER_SUFFIX=""
else
  ITER_SUFFIX="_n${N_ITER}"
fi
RESULT_TAG="w${WORKERS}${ITER_SUFFIX}_i${FAST_MAX_NMF_ITER}"
REF_ROOT="$DATA_ROOT/cnmf_optimized_w${WORKERS}${ITER_SUFFIX}"
FAST_PRE_ROOT="$DATA_ROOT/fastcnmf_independent_corecache_hvg_notp"
FAST_ROOT="$DATA_ROOT/fastcnmf_exact_cd_i${FAST_MAX_NMF_ITER}${ITER_SUFFIX}_corecache_hvg_notp_w${WORKERS}"
REF_RUN="${DATASET_ID}_cnmf_optimized_w${WORKERS}${ITER_SUFFIX}"
FAST_RUN="${DATASET_ID}_fastcnmf_exact_cd_i${FAST_MAX_NMF_ITER}${ITER_SUFFIX}_corecache_hvg_notp_w${WORKERS}"

REF_LOG="$REF_ROOT/logs"
FAST_PRE_LOG="$FAST_PRE_ROOT/logs"
FAST_LOG="$FAST_ROOT/logs"
mkdir -p "$REF_LOG" "$FAST_PRE_LOG" "$FAST_LOG"

stage_done() {
  local log="$1"
  [[ "$FORCE" == 0 && -f "$log" ]] && grep -q $'\tExit status: 0' "$log"
}

run_stage() {
  local name="$1"
  local log_dir="$2"
  shift 2
  local time_log="$log_dir/${name}.time.log"
  if stage_done "$time_log"; then
    echo "[skip] $name: $time_log"
    return 0
  fi
  echo "[run] $name"
  /usr/bin/time -v -o "$time_log" "$@" > "$log_dir/${name}.stdout.log" 2> "$log_dir/${name}.stderr.log"
}

echo "[dataset] $DATASET_ID"
echo "[input] $INPUT_H5AD"
echo "[out] $DATA_ROOT"

run_stage reference_preprocess "$REF_LOG" \
  "$REFERENCE_PYTHON" -m fastcnmf cnmf-preprocess \
    --input-h5ad "$INPUT_H5AD" \
    --output-prefix "$REF_ROOT/preprocess/cnmf_input" \
    --sample-key "$SAMPLE_KEY" \
    --n-top-genes "$N_TOP_GENES" \
    --theta 1 \
    --max-iter-harmony 20 \
    --seed "$SEED"

run_stage reference_plan_nmf "$REF_LOG" \
  bash -c "$(printf '%q ' "$CANDIDATE_PYTHON" -m fastcnmf plan-tasks --run-name "$REF_RUN" -k "${K_ARRAY[@]}" --n-iter "$N_ITER" --seed "$SEED" --output "$REF_ROOT/nmf/task_manifest.json") && $(printf '%q ' "$CNMF_BIN" prepare --output-dir "$REF_ROOT/nmf" --name "$REF_RUN" --counts "$REF_ROOT/preprocess/cnmf_input.Corrected.HVG.Varnorm.h5ad" --tpm "$REF_ROOT/preprocess/cnmf_input.TP10K.h5ad" --genes-file "$REF_ROOT/preprocess/cnmf_input.Corrected.HVGs.txt" -k "${K_ARRAY[@]}" --n-iter "$N_ITER" --seed "$SEED" --numgenes "$N_TOP_GENES" --max-nmf-iter "$REF_MAX_NMF_ITER")"

run_stage reference_factorize "$REF_LOG" \
  bash -c "set -euo pipefail; mkdir -p $(printf '%q' "$REF_ROOT/nmf/$REF_RUN/worker_logs"); status=0; for worker_i in \$(seq 0 $((WORKERS - 1))); do ($(printf '%q ' "$CNMF_BIN" factorize --output-dir "$REF_ROOT/nmf" --name "$REF_RUN" --total-workers "$WORKERS" --worker-index) \"\$worker_i\" > $(printf '%q' "$REF_ROOT/nmf/$REF_RUN/worker_logs")/factorize_worker_\${worker_i}.stdout.log 2> $(printf '%q' "$REF_ROOT/nmf/$REF_RUN/worker_logs")/factorize_worker_\${worker_i}.stderr.log) & done; for job in \$(jobs -p); do if ! wait \"\$job\"; then status=1; fi; done; exit \"\$status\""

run_stage reference_consensus "$REF_LOG" \
  bash -c "$(printf '%q ' "$CNMF_BIN" combine --output-dir "$REF_ROOT/nmf" --name "$REF_RUN") && $(printf '%q ' "$CNMF_BIN" consensus --output-dir "$REF_ROOT/nmf" --name "$REF_RUN") && $(printf '%q ' "$CNMF_BIN" k_selection_plot --output-dir "$REF_ROOT/nmf" --name "$REF_RUN")"

run_stage fast_preprocess "$FAST_PRE_LOG" \
  "$CANDIDATE_PYTHON" -m fastcnmf fast-preprocess \
    --input-h5ad "$INPUT_H5AD" \
    --output-prefix "$FAST_PRE_ROOT/preprocess/cnmf_input" \
    --sample-key "$SAMPLE_KEY" \
    --n-top-genes "$N_TOP_GENES" \
    --theta 1 \
    --max-iter-harmony 20 \
    --seed "$SEED" \
    --write-core-cache \
    --core-dtype float32 \
    --no-tp10k-h5ad

run_stage fast_plan_prepare "$FAST_LOG" \
  bash -c "$(printf '%q ' "$CANDIDATE_PYTHON" -m fastcnmf plan-tasks --run-name "$FAST_RUN" -k "${K_ARRAY[@]}" --n-iter "$N_ITER" --seed "$SEED" --output "$FAST_ROOT/nmf/task_manifest.json") && $(printf '%q ' "$CANDIDATE_PYTHON" -m fastcnmf fast-prepare --corrected-h5ad "$FAST_PRE_ROOT/preprocess/cnmf_input.Corrected.HVG.Varnorm.h5ad" --hvgs-txt "$FAST_PRE_ROOT/preprocess/cnmf_input.Corrected.HVGs.txt" --output-dir "$FAST_ROOT/nmf" --run-name "$FAST_RUN" -k "${K_ARRAY[@]}" --n-iter "$N_ITER" --seed "$SEED" --max-nmf-iter "$FAST_MAX_NMF_ITER" --norm-dtype float32 --norm-store npy --precomputed-norm-npy "$FAST_PRE_ROOT/preprocess/cnmf_input.NormCounts.float32.npy" --precomputed-norm-obs-names "$FAST_PRE_ROOT/preprocess/cnmf_input.NormCounts.obs.txt" --precomputed-norm-var-names "$FAST_PRE_ROOT/preprocess/cnmf_input.NormCounts.var.txt" --precomputed-tpm-stats "$FAST_PRE_ROOT/preprocess/cnmf_input.TP10K.stats.df.npz" --precomputed-tpm-hvg-raw "$FAST_PRE_ROOT/preprocess/cnmf_input.TP10K.HVG.raw.npz" --precomputed-tpm-hvg-scaled "$FAST_PRE_ROOT/preprocess/cnmf_input.TP10K.HVG.scaled.npz" --precomputed-tpm-hvg-obs-names "$FAST_PRE_ROOT/preprocess/cnmf_input.TP10K.HVG.obs.txt" --precomputed-tpm-hvg-var-names "$FAST_PRE_ROOT/preprocess/cnmf_input.TP10K.HVG.var.txt")"

run_stage fast_factorize "$FAST_LOG" \
  "$CANDIDATE_PYTHON" -m fastcnmf fast-factorize \
    --output-dir "$FAST_ROOT/nmf" \
    --run-name "$FAST_RUN" \
    --workers "$WORKERS" \
    --backend exact

run_stage fast_consensus_lite "$FAST_LOG" \
  "$CANDIDATE_PYTHON" -m fastcnmf fast-consensus \
    --output-dir "$FAST_ROOT/nmf" \
    --run-name "$FAST_RUN" \
    --workers "$CONSENSUS_WORKERS" \
    -k "${K_ARRAY[@]}" \
    --lite

run_stage compare_cnmf "$FAST_LOG" \
  "$CANDIDATE_PYTHON" -m fastcnmf compare-cnmf \
    --reference-dir "$REF_ROOT/nmf/$REF_RUN" \
    --candidate-dir "$FAST_ROOT/nmf/$FAST_RUN" \
    --reference-run-name "$REF_RUN" \
    --candidate-run-name "$FAST_RUN" \
    -k "${K_ARRAY[@]}" \
    --output-json "$DATA_ROOT/same_worker_${RESULT_TAG}_compare.json"

DATASET_ID="$DATASET_ID" DATA_ROOT="$DATA_ROOT" REF_ROOT="$REF_ROOT" FAST_PRE_ROOT="$FAST_PRE_ROOT" FAST_ROOT="$FAST_ROOT" WORKERS="$WORKERS" N_ITER="$N_ITER" RESULT_TAG="$RESULT_TAG" CANDIDATE_MAX_ITER="$FAST_MAX_NMF_ITER" "$CANDIDATE_PYTHON" - <<'PY'
import json
import os
import re
from pathlib import Path

dataset_id = os.environ["DATASET_ID"]
data_root = Path(os.environ["DATA_ROOT"])
ref_root = Path(os.environ["REF_ROOT"])
fast_pre_root = Path(os.environ["FAST_PRE_ROOT"])
fast_root = Path(os.environ["FAST_ROOT"])
workers = os.environ["WORKERS"]
n_iter = os.environ["N_ITER"]
result_tag = os.environ["RESULT_TAG"]
candidate_max_iter = os.environ["CANDIDATE_MAX_ITER"]

def parse_time(path: Path) -> dict:
    text = path.read_text(errors="replace")
    def grab(label: str):
        m = re.search(rf"{re.escape(label)}: (.+)", text)
        if not m:
            return None
        value = m.group(1).strip()
        if label == "Elapsed (wall clock) time (h:mm:ss or m:ss)":
            parts = value.split(":")
            if len(parts) == 3:
                return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
            if len(parts) == 2:
                return int(parts[0]) * 60 + float(parts[1])
            return float(parts[0])
        if label == "Maximum resident set size (kbytes)":
            return round(float(value) / 1024, 1)
        if label == "Percent of CPU this job got":
            return float(value.rstrip("%"))
        if label == "Exit status":
            return int(value)
        return float(value)
    return {
        "path": str(path),
        "elapsed_seconds": grab("Elapsed (wall clock) time (h:mm:ss or m:ss)"),
        "user_seconds": grab("User time (seconds)"),
        "system_seconds": grab("System time (seconds)"),
        "cpu_percent": grab("Percent of CPU this job got"),
        "max_rss_mb": grab("Maximum resident set size (kbytes)"),
        "exit_status": grab("Exit status"),
    }

def load_compare(path: Path) -> dict:
    d = json.loads(path.read_text())
    return {
        "path": str(path),
        "mean_spectra_cosine": d.get("mean_spectra_cosine"),
        "mean_usage_pearson": d.get("mean_usage_pearson"),
        "passes_95pct_gate": d.get("passes_95pct_gate"),
        "passes_all_k_95pct_gate": d.get("passes_all_k_95pct_gate"),
        "min_k_overall_consistency": d.get("min_k_overall_consistency"),
        "k_overall": {str(c["k"]): c.get("overall_consistency") for c in d.get("comparisons", [])},
    }

reference_steps = {
    "preprocess": parse_time(ref_root / "logs/reference_preprocess.time.log"),
    "prepare": parse_time(ref_root / "logs/reference_plan_nmf.time.log"),
    "factorize": parse_time(ref_root / "logs/reference_factorize.time.log"),
    "consensus": parse_time(ref_root / "logs/reference_consensus.time.log"),
}
fast_steps = {
    "preprocess": parse_time(fast_pre_root / "logs/fast_preprocess.time.log"),
    "prepare": parse_time(fast_root / "logs/fast_plan_prepare.time.log"),
    "factorize": parse_time(fast_root / "logs/fast_factorize.time.log"),
    "consensus": parse_time(fast_root / "logs/fast_consensus_lite.time.log"),
}
reference_total = sum(v["elapsed_seconds"] for v in reference_steps.values())
fast_total = sum(v["elapsed_seconds"] for v in fast_steps.values())
quality = load_compare(data_root / f"same_worker_{result_tag}_compare.json")
summary = {
    "dataset_id": dataset_id,
    "benchmark_date": "2026-06-15",
    "workers": int(workers),
    "n_iter": int(n_iter),
    "candidate_max_nmf_iter": int(candidate_max_iter),
    "reference_total_seconds": round(reference_total, 2),
    "fastcnmf_total_seconds": round(fast_total, 2),
    "speedup": round(reference_total / fast_total, 2),
    "reference_steps": reference_steps,
    "fastcnmf_steps": fast_steps,
    "quality": quality,
    "acceptance": {
        "target_speedup": 3.0,
        "min_consistency": 0.95,
        "passes_speed": reference_total / fast_total >= 3.0,
        "passes_quality": bool(quality.get("passes_all_k_95pct_gate")),
    },
}
json_path = data_root / f"same_worker_{result_tag}_fairness_summary.json"
json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
md_path = data_root / f"same_worker_{result_tag}_fairness_summary.md"
md_path.write_text(
    "\n".join(
        [
            f"# Same-worker benchmark summary: {dataset_id}",
            "",
            f"Workers: `{workers}`. NMF replicates per k: `{n_iter}`. Candidate max NMF iterations: `{candidate_max_iter}`.",
            "",
            "| lane | preprocess | prepare | factorize | consensus | total | speedup |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            f"| cNMF optimized | {reference_steps['preprocess']['elapsed_seconds']:.2f} | {reference_steps['prepare']['elapsed_seconds']:.2f} | {reference_steps['factorize']['elapsed_seconds']:.2f} | {reference_steps['consensus']['elapsed_seconds']:.2f} | {reference_total:.2f} | 1.00x |",
            f"| FastCNMF | {fast_steps['preprocess']['elapsed_seconds']:.2f} | {fast_steps['prepare']['elapsed_seconds']:.2f} | {fast_steps['factorize']['elapsed_seconds']:.2f} | {fast_steps['consensus']['elapsed_seconds']:.2f} | {fast_total:.2f} | {reference_total / fast_total:.2f}x |",
            "",
            f"Quality: mean spectra cosine `{quality['mean_spectra_cosine']:.4f}`, mean usage Pearson `{quality['mean_usage_pearson']:.4f}`, strict all-k gate `{quality['passes_all_k_95pct_gate']}`, min-k overall consistency `{quality['min_k_overall_consistency']:.4f}`.",
            "",
        ]
    ),
    encoding="utf-8",
)
latest_json = data_root / f"same_worker_w{workers}_fairness_summary.json"
latest_md = data_root / f"same_worker_w{workers}_fairness_summary.md"
latest_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
latest_md.write_text(md_path.read_text(encoding="utf-8"), encoding="utf-8")
print(json.dumps({"summary_json": str(json_path), "summary_md": str(md_path), "speedup": round(reference_total / fast_total, 2), "quality_pass": quality.get("passes_all_k_95pct_gate")}, indent=2))
PY
