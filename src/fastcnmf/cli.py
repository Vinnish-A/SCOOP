from __future__ import annotations

import argparse
import json
from pathlib import Path

from .benchmark_manifest import FastCNMFBenchmarkManifest, build_default_manifest
from .cnmf_reference import run_cnmf_reference_preprocess
from .fast_consensus import run_fast_consensus
from .fast_factorize import run_fast_factorize
from .fast_prepare import run_fast_prepare
from .preprocess import run_fast_harmony_preprocess
from .preprocess_quality import compare_preprocess_outputs
from .quality import compare_cnmf_dirs, evaluate_benchmark_gate
from .resources import choose_replicate_batch_size, estimate_nmf_memory
from .run_bundle import write_run_bundle
from .runner import build_execution_plan
from .stage_report import write_stage_report
from .tasks import build_nmf_tasks


def cmd_status(args: argparse.Namespace) -> int:
    report = Path(args.report)
    if not report.exists():
        raise FileNotFoundError(f"status report not found: {report}")
    print(report.read_text(encoding="utf-8"))
    return 0


def cmd_plan_tasks(args: argparse.Namespace) -> int:
    manifest = build_nmf_tasks(
        run_name=args.run_name,
        k_values=tuple(args.k),
        n_iter=args.n_iter,
        seed=args.seed,
    )
    manifest.to_json(Path(args.output))
    print(json.dumps({"tasks": len(manifest.tasks), "output": args.output}, indent=2))
    return 0


def cmd_estimate_memory(args: argparse.Namespace) -> int:
    estimate = estimate_nmf_memory(
        observations=args.observations,
        genes=args.genes,
        components=args.components,
        dtype=args.dtype,
        replicate_batch_size=args.replicate_batch_size,
    )
    batch_size = choose_replicate_batch_size(
        observations=args.observations,
        genes=args.genes,
        components=args.components,
        available_bytes=int(args.available_gib * 1024**3),
        dtype=args.dtype,
        safety_fraction=args.safety_fraction,
        max_batch_size=args.max_batch_size,
    )
    print(
        json.dumps(
            {
                "estimate_gib": round(estimate.estimated_gib, 3),
                "estimate_bytes": estimate.estimated_bytes,
                "recommended_replicate_batch_size": batch_size,
            },
            indent=2,
        )
    )
    return 0


def cmd_benchmark_gate(args: argparse.Namespace) -> int:
    result = evaluate_benchmark_gate(
        summary_json=Path(args.summary_json),
        reference_dir=Path(args.reference_dir),
        candidate_dir=Path(args.candidate_dir),
        run_name=args.run_name,
        k_values=tuple(args.k),
        min_time_saved_fraction=args.min_time_saved_fraction,
        min_consistency=args.min_consistency,
    )
    output_json = Path(args.output_json)
    result.to_json(output_json)
    if args.output_md:
        write_gate_markdown(result=result, path=Path(args.output_md))
    print(json.dumps({"passed": result.passed, "output_json": str(output_json)}, indent=2))
    return 0 if result.passed else 1


def cmd_benchmark_manifest(args: argparse.Namespace) -> int:
    manifest = build_default_manifest(root=Path(args.root), output_root=Path(args.output_root))
    output = Path(args.output)
    manifest.to_json(output)
    print(
        json.dumps(
            {
                "output": str(output),
                "datasets": len(manifest.datasets),
                "missing_tiers": list(manifest.missing_tiers),
                "target_speedup": manifest.fairness.target_speedup,
            },
            indent=2,
        )
    )
    return 0


def cmd_plan_run(args: argparse.Namespace) -> int:
    manifest = FastCNMFBenchmarkManifest.from_json(Path(args.manifest))
    plan = build_execution_plan(manifest)
    output = Path(args.output)
    plan.to_json(output)
    print(json.dumps({"output": str(output), "stages": len(plan.stages)}, indent=2))
    return 0


def cmd_write_run_bundle(args: argparse.Namespace) -> int:
    manifest = FastCNMFBenchmarkManifest.from_json(Path(args.manifest))
    plan = build_execution_plan(manifest)
    bundle = write_run_bundle(
        manifest=manifest,
        plan=plan,
        output_dir=Path(args.output_dir),
        reference_python=args.reference_python,
        candidate_python=args.candidate_python,
        reference_cuda_visible_devices=args.reference_cuda_visible_devices,
        candidate_cuda_visible_devices=args.candidate_cuda_visible_devices,
        profile=args.profile,
        blas_threads=args.blas_threads,
    )
    output = Path(args.output_dir) / "run_bundle.json"
    bundle.to_json(output)
    print(
        json.dumps(
            {
                "output": str(output),
                "scripts": len(bundle.scripts),
                "implemented_scripts": sum(1 for script in bundle.scripts if script.implemented),
            },
            indent=2,
        )
    )
    return 0


def cmd_cnmf_preprocess(args: argparse.Namespace) -> int:
    manifest = run_cnmf_reference_preprocess(
        input_h5ad=Path(args.input_h5ad),
        output_prefix=Path(args.output_prefix),
        sample_key=args.sample_key,
        n_top_genes=args.n_top_genes,
        librarysize_targetsum=args.librarysize_targetsum,
        theta=args.theta,
        max_iter_harmony=args.max_iter_harmony,
        seed=args.seed,
    )
    print(json.dumps(manifest, indent=2))
    return 0


def cmd_fast_preprocess(args: argparse.Namespace) -> int:
    manifest = run_fast_harmony_preprocess(
        input_h5ad=Path(args.input_h5ad),
        output_prefix=Path(args.output_prefix),
        sample_key=args.sample_key,
        n_top_genes=args.n_top_genes,
        librarysize_targetsum=args.librarysize_targetsum,
        theta=args.theta,
        lamb=args.lamb,
        max_iter_harmony=args.max_iter_harmony,
        seed=args.seed,
        write_core_cache=args.write_core_cache,
        core_dtype=args.core_dtype,
        write_corrected_h5ad=not args.no_corrected_h5ad,
        write_tp10k_h5ad=not args.no_tp10k_h5ad,
    )
    print(json.dumps(manifest, indent=2))
    return 0


def cmd_compare_preprocess(args: argparse.Namespace) -> int:
    result = compare_preprocess_outputs(
        reference_h5ad=Path(args.reference_h5ad),
        candidate_h5ad=Path(args.candidate_h5ad),
        output_json=Path(args.output_json),
        chunk_size=args.chunk_size,
    )
    print(
        json.dumps(
            {
                "output_json": args.output_json,
                "cosine": result["matrix"]["cosine"],
                "pearson": result["matrix"]["pearson"],
                "passes_95pct_input_gate": result["passes_95pct_input_gate"],
            },
            indent=2,
        )
    )
    return 0 if result["passes_95pct_input_gate"] else 1


def cmd_stage_report(args: argparse.Namespace) -> int:
    report = write_stage_report(
        dataset_id=args.dataset_id,
        lane_id=args.lane_id,
        lane_root=Path(args.lane_root),
        logs_dir=Path(args.logs_dir),
        output_json=Path(args.output_json),
    )
    print(json.dumps({"output_json": args.output_json, "stages": len(report["stages"])}, indent=2))
    return 0


def cmd_compare_cnmf(args: argparse.Namespace) -> int:
    result = compare_cnmf_dirs(
        reference_dir=Path(args.reference_dir),
        candidate_dir=Path(args.candidate_dir),
        reference_run_name=args.reference_run_name,
        candidate_run_name=args.candidate_run_name,
        k_values=tuple(args.k),
        output_json=Path(args.output_json),
    )
    print(
        json.dumps(
            {
                "output_json": args.output_json,
                "mean_spectra_cosine": result["mean_spectra_cosine"],
                "mean_usage_pearson": result["mean_usage_pearson"],
                "passes_95pct_gate": result["passes_95pct_gate"],
            },
            indent=2,
        )
    )
    return 0 if result["passes_95pct_gate"] else 1


def cmd_fast_factorize(args: argparse.Namespace) -> int:
    result = run_fast_factorize(
        output_dir=Path(args.output_dir),
        run_name=args.run_name,
        workers=args.workers,
        skip_completed=args.skip_completed_runs,
        compressed=args.compressed,
        backend=args.backend,
        minibatch_batch_size=args.minibatch_batch_size,
        mu_dtype=args.mu_dtype,
    )
    print(json.dumps(result, indent=2))
    return 0


def cmd_fast_prepare(args: argparse.Namespace) -> int:
    result = run_fast_prepare(
        corrected_h5ad=Path(args.corrected_h5ad),
        tpm_h5ad=Path(args.tpm_h5ad) if args.tpm_h5ad else None,
        hvgs_txt=Path(args.hvgs_txt),
        output_dir=Path(args.output_dir),
        run_name=args.run_name,
        k_values=tuple(args.k),
        n_iter=args.n_iter,
        seed=args.seed,
        beta_loss=args.beta_loss,
        max_nmf_iter=args.max_nmf_iter,
        init=args.init,
        norm_dtype=args.norm_dtype,
        norm_store=args.norm_store,
        precomputed_norm_npy=Path(args.precomputed_norm_npy) if args.precomputed_norm_npy else None,
        precomputed_norm_obs_names=Path(args.precomputed_norm_obs_names) if args.precomputed_norm_obs_names else None,
        precomputed_norm_var_names=Path(args.precomputed_norm_var_names) if args.precomputed_norm_var_names else None,
        precomputed_tpm_stats=Path(args.precomputed_tpm_stats) if args.precomputed_tpm_stats else None,
        precomputed_tpm_hvg_raw=Path(args.precomputed_tpm_hvg_raw) if args.precomputed_tpm_hvg_raw else None,
        precomputed_tpm_hvg_scaled=Path(args.precomputed_tpm_hvg_scaled) if args.precomputed_tpm_hvg_scaled else None,
        precomputed_tpm_hvg_obs_names=Path(args.precomputed_tpm_hvg_obs_names) if args.precomputed_tpm_hvg_obs_names else None,
        precomputed_tpm_hvg_var_names=Path(args.precomputed_tpm_hvg_var_names) if args.precomputed_tpm_hvg_var_names else None,
    )
    print(json.dumps(result, indent=2))
    return 0


def cmd_fast_consensus(args: argparse.Namespace) -> int:
    result = run_fast_consensus(
        output_dir=Path(args.output_dir),
        run_name=args.run_name,
        workers=args.workers,
        k_values=tuple(args.k) if args.k else None,
        density_threshold=args.local_density_threshold,
        local_neighborhood_size=args.local_neighborhood_size,
        show_clustering=args.show_clustering,
        build_reference=not args.no_build_reference,
        k_selection=not args.no_k_selection,
        lite=args.lite,
    )
    print(json.dumps(result, indent=2))
    return 0


def write_gate_markdown(result, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for k_result in result.consistency:
        rows.append(
            "| {k} | {overall:.6f} | {spectra:.6f} | {usage:.6f} | {jaccard:.6f} |".format(
                k=k_result.k,
                overall=k_result.overall_consistency,
                spectra=k_result.mean_spectra_cosine,
                usage=k_result.mean_usage_pearson,
                jaccard=k_result.mean_top_gene_jaccard,
            )
        )
    path.write_text(
        "\n".join(
            [
                "# FastCNMF Benchmark Gate",
                "",
                f"Passed: `{result.passed}`",
                "",
                "## Speed",
                "",
                f"- Reference factorize seconds: `{result.reference_factorize_seconds}`",
                f"- Candidate factorize seconds: `{result.candidate_factorize_seconds}`",
                f"- Speedup: `{result.speedup:.3f}x`",
                f"- Time saved fraction: `{result.time_saved_fraction:.3%}`",
                f"- Required time saved fraction: `{result.min_required_time_saved_fraction:.3%}`",
                "",
                "## Consistency",
                "",
                f"Required minimum consistency: `{result.min_required_consistency:.3%}`",
                "",
                "| k | overall_consistency | mean_spectra_cosine | mean_usage_pearson | mean_top_gene_jaccard |",
                "| ---: | ---: | ---: | ---: | ---: |",
                *rows,
                "",
                "Overall consistency is the stricter of mean matched spectra cosine and mean matched usage Pearson.",
                "",
                "## JSON",
                "",
                "The full program matching details are available in the sibling JSON report.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="fastcnmf")
    sub = parser.add_subparsers(dest="command", required=True)

    status = sub.add_parser("status", help="print the FastCNMF execution status document")
    status.add_argument("--report", default="docs/fastcnmf/02_execution_status.md")
    status.set_defaults(func=cmd_status)

    plan_tasks = sub.add_parser("plan-tasks", help="write an NMF replicate task manifest")
    plan_tasks.add_argument("--run-name", default="fastcnmf")
    plan_tasks.add_argument("-k", type=int, nargs="+", required=True)
    plan_tasks.add_argument("--n-iter", type=int, required=True)
    plan_tasks.add_argument("--seed", type=int, default=20260614)
    plan_tasks.add_argument("--output", required=True)
    plan_tasks.set_defaults(func=cmd_plan_tasks)

    memory = sub.add_parser("estimate-memory", help="estimate dense NMF memory usage")
    memory.add_argument("--observations", type=int, required=True)
    memory.add_argument("--genes", type=int, required=True)
    memory.add_argument("--components", type=int, required=True)
    memory.add_argument("--dtype", default="float32")
    memory.add_argument("--replicate-batch-size", type=int, default=1)
    memory.add_argument("--available-gib", type=float, default=24)
    memory.add_argument("--safety-fraction", type=float, default=0.8)
    memory.add_argument("--max-batch-size", type=int, default=8)
    memory.set_defaults(func=cmd_estimate_memory)

    gate = sub.add_parser("benchmark-gate", help="check speed and consistency benchmark gates")
    gate.add_argument("--summary-json", required=True)
    gate.add_argument("--reference-dir", required=True)
    gate.add_argument("--candidate-dir", required=True)
    gate.add_argument("--run-name", required=True)
    gate.add_argument("-k", type=int, nargs="+", required=True)
    gate.add_argument("--min-time-saved-fraction", type=float, default=0.50)
    gate.add_argument("--min-consistency", type=float, default=0.95)
    gate.add_argument("--output-json", required=True)
    gate.add_argument("--output-md", default=None)
    gate.set_defaults(func=cmd_benchmark_gate)

    manifest = sub.add_parser("benchmark-manifest", help="write the S1/S2 cold-start benchmark manifest")
    manifest.add_argument("--root", default=".")
    manifest.add_argument("--output-root", default="tmp/fastcnmf_large_benchmark")
    manifest.add_argument("--output", required=True)
    manifest.set_defaults(func=cmd_benchmark_manifest)

    plan_run = sub.add_parser("plan-run", help="expand a benchmark manifest into independent execution stages")
    plan_run.add_argument("--manifest", required=True)
    plan_run.add_argument("--output", required=True)
    plan_run.set_defaults(func=cmd_plan_run)

    bundle = sub.add_parser("write-run-bundle", help="write executable shell scripts for a benchmark manifest")
    bundle.add_argument("--manifest", required=True)
    bundle.add_argument("--output-dir", required=True)
    bundle.add_argument("--reference-python", default="python")
    bundle.add_argument("--candidate-python", default="python")
    bundle.add_argument("--reference-cuda-visible-devices", default=None)
    bundle.add_argument("--candidate-cuda-visible-devices", default=None)
    bundle.add_argument("--profile", choices=["smoke", "production"], default="production")
    bundle.add_argument("--blas-threads", type=int, default=1)
    bundle.set_defaults(func=cmd_write_run_bundle)

    cnmf_pre = sub.add_parser("cnmf-preprocess", help="run cNMF reference preprocessing from a cold H5AD")
    cnmf_pre.add_argument("--input-h5ad", required=True)
    cnmf_pre.add_argument("--output-prefix", required=True)
    cnmf_pre.add_argument("--sample-key", default="sample_id")
    cnmf_pre.add_argument("--n-top-genes", type=int, default=3000)
    cnmf_pre.add_argument("--librarysize-targetsum", type=float, default=1e4)
    cnmf_pre.add_argument("--theta", type=float, default=1.0)
    cnmf_pre.add_argument("--max-iter-harmony", type=int, default=20)
    cnmf_pre.add_argument("--seed", type=int, default=20260614)
    cnmf_pre.set_defaults(func=cmd_cnmf_preprocess)

    fast_pre = sub.add_parser("fast-preprocess", help="run FastCNMF independent preprocessing from a cold H5AD")
    fast_pre.add_argument("--input-h5ad", required=True)
    fast_pre.add_argument("--output-prefix", required=True)
    fast_pre.add_argument("--sample-key", default="sample_id")
    fast_pre.add_argument("--n-top-genes", type=int, default=3000)
    fast_pre.add_argument("--librarysize-targetsum", type=float, default=1e4)
    fast_pre.add_argument("--theta", type=float, default=1.0)
    fast_pre.add_argument("--lamb", type=float, default=1.0)
    fast_pre.add_argument("--max-iter-harmony", type=int, default=20)
    fast_pre.add_argument("--seed", type=int, default=20260614)
    fast_pre.add_argument("--write-core-cache", action="store_true", default=False)
    fast_pre.add_argument("--core-dtype", choices=["float64", "float32"], default="float32")
    fast_pre.add_argument("--no-corrected-h5ad", action="store_true", default=False)
    fast_pre.add_argument("--no-tp10k-h5ad", action="store_true", default=False)
    fast_pre.set_defaults(func=cmd_fast_preprocess)

    cmp_pre = sub.add_parser("compare-preprocess", help="compare reference and candidate corrected preprocess H5AD files")
    cmp_pre.add_argument("--reference-h5ad", required=True)
    cmp_pre.add_argument("--candidate-h5ad", required=True)
    cmp_pre.add_argument("--output-json", required=True)
    cmp_pre.add_argument("--chunk-size", type=int, default=2000)
    cmp_pre.set_defaults(func=cmd_compare_preprocess)

    stage_report = sub.add_parser("stage-report", help="summarize a lane's stage time logs and output size")
    stage_report.add_argument("--dataset-id", required=True)
    stage_report.add_argument("--lane-id", required=True)
    stage_report.add_argument("--lane-root", required=True)
    stage_report.add_argument("--logs-dir", required=True)
    stage_report.add_argument("--output-json", required=True)
    stage_report.set_defaults(func=cmd_stage_report)

    cmp_cnmf = sub.add_parser("compare-cnmf", help="compare two cNMF output directories with different run names")
    cmp_cnmf.add_argument("--reference-dir", required=True)
    cmp_cnmf.add_argument("--candidate-dir", required=True)
    cmp_cnmf.add_argument("--reference-run-name", required=True)
    cmp_cnmf.add_argument("--candidate-run-name", required=True)
    cmp_cnmf.add_argument("-k", type=int, nargs="+", required=True)
    cmp_cnmf.add_argument("--output-json", required=True)
    cmp_cnmf.set_defaults(func=cmd_compare_cnmf)

    fast_prepare = sub.add_parser("fast-prepare", help="write cNMF-compatible prepare artifacts directly")
    fast_prepare.add_argument("--corrected-h5ad", required=True)
    fast_prepare.add_argument("--tpm-h5ad", default=None)
    fast_prepare.add_argument("--hvgs-txt", required=True)
    fast_prepare.add_argument("--output-dir", required=True)
    fast_prepare.add_argument("--run-name", required=True)
    fast_prepare.add_argument("-k", type=int, nargs="+", required=True)
    fast_prepare.add_argument("--n-iter", type=int, required=True)
    fast_prepare.add_argument("--seed", type=int, default=20260614)
    fast_prepare.add_argument("--beta-loss", choices=["frobenius", "kullback-leibler", "itakura-saito"], default="frobenius")
    fast_prepare.add_argument("--max-nmf-iter", type=int, default=50)
    fast_prepare.add_argument("--init", choices=["random", "nndsvd"], default="random")
    fast_prepare.add_argument("--norm-dtype", choices=["float64", "float32"], default="float64")
    fast_prepare.add_argument("--norm-store", choices=["h5ad", "npy", "both"], default="h5ad")
    fast_prepare.add_argument("--precomputed-norm-npy", default=None)
    fast_prepare.add_argument("--precomputed-norm-obs-names", default=None)
    fast_prepare.add_argument("--precomputed-norm-var-names", default=None)
    fast_prepare.add_argument("--precomputed-tpm-stats", default=None)
    fast_prepare.add_argument("--precomputed-tpm-hvg-raw", default=None)
    fast_prepare.add_argument("--precomputed-tpm-hvg-scaled", default=None)
    fast_prepare.add_argument("--precomputed-tpm-hvg-obs-names", default=None)
    fast_prepare.add_argument("--precomputed-tpm-hvg-var-names", default=None)
    fast_prepare.set_defaults(func=cmd_fast_prepare)

    fast_factorize = sub.add_parser("fast-factorize", help="run FastCNMF scheduled cNMF-compatible NMF replicates")
    fast_factorize.add_argument("--output-dir", required=True)
    fast_factorize.add_argument("--run-name", required=True)
    fast_factorize.add_argument("--workers", type=int, default=4)
    fast_factorize.add_argument("--skip-completed-runs", action="store_true", default=False)
    fast_factorize.add_argument("--compressed", action="store_true", default=False)
    fast_factorize.add_argument("--backend", choices=["exact", "minibatch", "cupy-mu"], default="exact")
    fast_factorize.add_argument("--minibatch-batch-size", type=int, default=1024)
    fast_factorize.add_argument("--mu-dtype", choices=["float32", "float64"], default="float32")
    fast_factorize.set_defaults(func=cmd_fast_factorize)

    fast_consensus = sub.add_parser("fast-consensus", help="run cNMF-compatible consensus with FastCNMF scheduling")
    fast_consensus.add_argument("--output-dir", required=True)
    fast_consensus.add_argument("--run-name", required=True)
    fast_consensus.add_argument("--workers", type=int, default=4)
    fast_consensus.add_argument("-k", type=int, nargs="+", default=None)
    fast_consensus.add_argument("--local-density-threshold", type=float, default=0.5)
    fast_consensus.add_argument("--local-neighborhood-size", type=float, default=0.30)
    fast_consensus.add_argument("--show-clustering", action="store_true", default=False)
    fast_consensus.add_argument("--no-build-reference", action="store_true", default=False)
    fast_consensus.add_argument("--no-k-selection", action="store_true", default=False)
    fast_consensus.add_argument("--lite", action="store_true", default=False)
    fast_consensus.set_defaults(func=cmd_fast_consensus)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
