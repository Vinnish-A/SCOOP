from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from fastde.abundance import run_abundance
from fastde.abundance_metrics import concordance_index


def make_fixture(output_dir: Path, n_samples: int, cells_per_sample: int, n_celltypes: int, seed: int) -> dict[str, str]:
    rng = np.random.default_rng(seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    celltypes = ["RiskCells"] + [f"Other{i}" for i in range(1, n_celltypes)]
    samples = [f"S{i:03d}" for i in range(n_samples)]
    obs_rows = []
    x_rows = []
    truth_risk = []
    counts = []
    for i, sample in enumerate(samples):
        risk_prop = 0.08 + 0.72 * (i / max(n_samples - 1, 1))
        truth_risk.append(risk_prop)
        alpha = np.ones(n_celltypes)
        alpha[0] = 1 + risk_prop * 12
        probs = rng.dirichlet(alpha)
        probs[0] = risk_prop
        probs[1:] = probs[1:] / probs[1:].sum() * (1.0 - risk_prop)
        sample_counts = rng.multinomial(cells_per_sample, probs)
        counts.append(sample_counts)
        for celltype, n_cells in zip(celltypes, sample_counts):
            for _ in range(int(n_cells)):
                feature = np.zeros(n_celltypes, dtype=float)
                feature[celltypes.index(celltype)] = 1.0
                feature = np.clip(feature + rng.normal(0, 0.02, size=n_celltypes), 0.0, None)
                x_rows.append(feature)
                obs_rows.append({"sample_id": sample, "cell_type_lvl3": celltype})

    obs = pd.DataFrame(obs_rows)
    x = np.asarray(x_rows, dtype=np.float32)
    adata = ad.AnnData(X=x.copy(), obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_celltypes)]))
    adata.obsm["X_abund_feat"] = x.copy()
    h5ad_path = output_dir / "synthetic_scsurvival_abundance.h5ad"
    adata.write_h5ad(h5ad_path)

    truth = np.asarray(truth_risk)
    time_values = np.maximum(1.0, 120.0 - 90.0 * truth + rng.normal(0.0, 3.0, size=n_samples))
    event_values = np.ones(n_samples, dtype=int)
    event_values[rng.choice(n_samples, size=max(1, n_samples // 5), replace=False)] = 0
    metadata = pd.DataFrame({"sample_id": samples, "OS_time": time_values, "OS_event": event_values, "truth_risk": truth})
    metadata_path = output_dir / "sample_metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    pd.DataFrame(counts, index=samples, columns=celltypes).to_csv(output_dir / "sample_by_celltype_counts.tsv", sep="\t")
    return {"h5ad": str(h5ad_path), "metadata": str(metadata_path)}


def peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return float(rss / 1024.0)


def run_fastde_child(args: argparse.Namespace) -> None:
    start = time.perf_counter()
    result = run_abundance(
        mode="survival",
        input_h5ad=args.input_h5ad,
        metadata=args.metadata,
        sample_key="sample_id",
        celltype_key="cell_type_lvl3",
        time_col="OS_time",
        event_col="OS_event",
        max_epochs=args.epochs,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )
    predictions = result.predictions.rename(columns={"risk": "fastde_risk"})
    top_cell = str(result.results.sort_values("pvalue", na_position="last").iloc[0]["cell_type"])
    metrics = {
        "framework": "FastDE abundance",
        "wall_seconds": time.perf_counter() - start,
        "peak_rss_mb": peak_rss_mb(),
        "cindex": result.metrics.get("concordance_index"),
        "top_cell_type": top_cell,
        "n_samples": result.manifest["n_samples"],
        "n_celltypes": result.manifest["n_celltypes"],
    }
    predictions.to_csv(Path(args.output_dir) / "fastde_reference_predictions.tsv", sep="\t", index=False)
    Path(args.output_dir, "fastde_reference_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")


def run_scsurvival_child(args: argparse.Namespace) -> None:
    start = time.perf_counter()
    from scSurvival.scsurvival import scSurvivalRun

    adata = ad.read_h5ad(args.input_h5ad)
    metadata = pd.read_csv(args.metadata, sep="\t").set_index("sample_id")
    surv = metadata[["OS_time", "OS_event"]].rename(columns={"OS_time": "time", "OS_event": "status"})
    adata_out, surv_out, _ = scSurvivalRun(
        adata,
        sample_column="sample_id",
        surv=surv,
        feature_flavor="Custom",
        feature_key="X_abund_feat",
        num_heads=None,
        dropout=0.0,
        predict_nMC=1,
        epochs=args.epochs,
        pretrain_epochs=0,
        lr=args.learning_rate,
        instance_batch_size=args.instance_batch_size,
        pretrain_batch_size=args.instance_batch_size,
        lambdas=(0.0, 0.0),
        entropy_threshold=1.0,
        patience=max(10, args.epochs),
        once_load_to_gpu=False,
        use_amp=False,
    )
    pred = surv_out[["time", "status", "patient_hazards"]].copy()
    pred["sample_id"] = pred.index.astype(str)
    cindex = concordance_index(pred["time"].to_numpy(), pred["patient_hazards"].to_numpy(), pred["status"].to_numpy())
    cell_scores = adata_out.obs.groupby("cell_type_lvl3", observed=True)["hazard_adj"].mean().sort_values(ascending=False)
    metrics = {
        "framework": "scSurvival reference",
        "wall_seconds": time.perf_counter() - start,
        "peak_rss_mb": peak_rss_mb(),
        "cindex": cindex,
        "top_cell_type": str(cell_scores.index[0]),
        "n_samples": int(pred.shape[0]),
        "n_cells": int(adata_out.n_obs),
    }
    out = Path(args.output_dir)
    pred.rename(columns={"patient_hazards": "scsurvival_risk"}).to_csv(out / "scsurvival_reference_predictions.tsv", sep="\t", index=False)
    cell_scores.rename("mean_hazard_adj").reset_index().to_csv(out / "scsurvival_celltype_scores.tsv", sep="\t", index=False)
    Path(out, "scsurvival_reference_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")


def run_child(python: str, script: Path, runner: str, args: argparse.Namespace, output_dir: Path, env_extra: dict[str, str] | None = None) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{Path.cwd() / 'src'}{os.pathsep}{env.get('PYTHONPATH', '')}"
    if env_extra:
        env.update(env_extra)
    if args.force_cpu:
        env["CUDA_VISIBLE_DEVICES"] = ""
    cmd = [
        python,
        str(script),
        "--runner",
        runner,
        "--input-h5ad",
        str(output_dir / "synthetic_scsurvival_abundance.h5ad"),
        "--metadata",
        str(output_dir / "sample_metadata.tsv"),
        "--output-dir",
        str(output_dir / runner),
        "--epochs",
        str(args.epochs),
        "--learning-rate",
        str(args.learning_rate),
        "--instance-batch-size",
        str(args.instance_batch_size),
    ]
    (output_dir / runner).mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True, env=env)


def summarize(output_dir: Path) -> None:
    fast = json.loads((output_dir / "fastde" / "fastde_reference_metrics.json").read_text(encoding="utf-8"))
    scs_path = output_dir / "scsurvival" / "scsurvival_reference_metrics.json"
    rows = [fast]
    consistency: dict[str, float | str | None] = {}
    if scs_path.exists():
        scs = json.loads(scs_path.read_text(encoding="utf-8"))
        rows.append(scs)
        f_pred = pd.read_csv(output_dir / "fastde" / "fastde_reference_predictions.tsv", sep="\t")
        s_pred = pd.read_csv(output_dir / "scsurvival" / "scsurvival_reference_predictions.tsv", sep="\t")
        merged = f_pred.merge(s_pred[["sample_id", "scsurvival_risk"]], on="sample_id")
        corr, pvalue = stats.spearmanr(merged["fastde_risk"], merged["scsurvival_risk"])
        consistency = {
            "risk_spearman": float(corr),
            "risk_spearman_pvalue": float(pvalue),
            "top_cell_type_match": str(fast.get("top_cell_type")) == str(scs.get("top_cell_type")),
            "fastde_top_cell_type": str(fast.get("top_cell_type")),
            "scsurvival_top_cell_type": str(scs.get("top_cell_type")),
            "speedup_vs_scsurvival": float(scs["wall_seconds"] / fast["wall_seconds"]) if fast["wall_seconds"] else None,
            "rss_ratio_scsurvival_vs_fastde": float(scs["peak_rss_mb"] / fast["peak_rss_mb"]) if fast["peak_rss_mb"] else None,
        }
    pd.DataFrame(rows).to_csv(output_dir / "benchmark_summary.tsv", sep="\t", index=False)
    (output_dir / "benchmark_consistency.json").write_text(json.dumps(consistency, indent=2), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", choices=["all", "fastde", "scsurvival"], default="all")
    parser.add_argument("--output-dir", default="tmp/fastde_abundance_scsurvival_reference")
    parser.add_argument("--input-h5ad", default=None)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--n-samples", type=int, default=24)
    parser.add_argument("--cells-per-sample", type=int, default=120)
    parser.add_argument("--n-celltypes", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--learning-rate", type=float, default=0.02)
    parser.add_argument("--instance-batch-size", type=int, default=256)
    parser.add_argument("--fastde-python", default=sys.executable)
    parser.add_argument("--scsurvival-python", default=None)
    parser.add_argument("--scsurvival-source", default="/tmp/scSurvival_ref")
    parser.add_argument("--force-cpu", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out = Path(args.output_dir)
    if args.runner == "fastde":
        run_fastde_child(args)
        return 0
    if args.runner == "scsurvival":
        run_scsurvival_child(args)
        return 0

    make_fixture(out, args.n_samples, args.cells_per_sample, args.n_celltypes, args.seed)
    script = Path(__file__).resolve()
    run_child(args.fastde_python, script, "fastde", args, out)
    if args.scsurvival_python:
        run_child(
            args.scsurvival_python,
            script,
            "scsurvival",
            args,
            out,
            env_extra={"PYTHONPATH": f"{Path.cwd() / 'src'}{os.pathsep}{args.scsurvival_source}"},
        )
    summarize(out)
    print(f"wrote benchmark summary to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
