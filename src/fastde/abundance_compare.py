"""Donor-level composition tests and held-out binary attention MIL.

MIL is a scSurvival-style composition-only adaptation: one-hot cell identities,
gated attention and a trained bag-level binary head. It is NOT expression-based
scSurvival, survival inference, or a cell-type annotation model.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import hashlib
import json
import math
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
import torch
from torch import nn

from .abundance_mil import GatedAttention


class BinaryCountBagMIL(nn.Module):
    """Exact compressed evaluation of repeated one-hot cell instances."""

    def __init__(self, n_types: int, hidden: int = 16):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(n_types, hidden), nn.ReLU(),
                                     nn.Linear(hidden, hidden), nn.ReLU())
        self.attention = GatedAttention(hidden)
        self.head = nn.Sequential(nn.Linear(hidden, hidden // 2), nn.ReLU(),
                                  nn.Dropout(0.1), nn.Linear(hidden // 2, 1))
        self.register_buffer("eye", torch.eye(n_types))

    def forward(self, counts):
        h = self.encoder(self.eye)
        # log multiplicity exactly reproduces cell-level softmax attention;
        # no rare cells are dropped and bag size alone supplies no signal.
        weights = torch.softmax(counts.log() + self.attention(h).flatten(), dim=1)
        return self.head(weights @ h).flatten()


def validate_inputs(counts: pd.DataFrame, metadata: pd.DataFrame):
    if not counts.index.is_unique or not metadata.index.is_unique:
        raise ValueError("Sample identifiers must be unique")
    if set(counts.index) != set(metadata.index):
        raise ValueError("Counts and metadata sample identities must match exactly")
    metadata = metadata.loc[counts.index]
    if metadata.donor_id.isna().any() or not metadata.donor_id.is_unique:
        raise ValueError("Use one bag per donor within each tissue; duplicate/missing donors rejected")
    if metadata.ACPA.isna().any() or set(metadata.ACPA) != {"Positive", "Negative"}:
        raise ValueError("Only ACPA Positive and Negative RA are allowed")
    x = counts.to_numpy(dtype=float)
    if not np.isfinite(x).all() or (x < 0).any() or not np.equal(x, np.floor(x)).all():
        raise ValueError("Counts must be finite, nonnegative integers")
    if (x.sum(axis=1) == 0).any():
        raise ValueError("Empty donor bags cannot be compared")
    return x, (metadata.ACPA == "Positive").to_numpy(dtype=int), metadata


def donor_tests(x, y, target, bootstrap=10000, seed=0):
    rng = np.random.default_rng(seed)
    target_n, totals = x[:, target], x.sum(axis=1)
    proportions = target_n / totals
    pos, neg = y == 1, y == 0
    a, b = proportions[pos], proportions[neg]
    delta = a.mean() - b.mean()
    pooled = target_n.sum() / totals.sum()
    roe_pos = (target_n[pos].sum() / totals[pos].sum()) / pooled
    roe_neg = (target_n[neg].sum() / totals[neg].sum()) / pooled
    roe_delta = roe_pos - roe_neg
    n, npos = len(y), int(pos.sum())
    # Exact donor-label permutations when tractable; never cell-level labels.
    exact = math.comb(n, npos) <= 200000
    assignments = (combinations(range(n), npos) if exact else
                   (rng.choice(n, npos, replace=False) for _ in range(100000)))
    p_delta = p_roe = permutations = 0
    for ids in assignments:
        ids = np.asarray(ids, dtype=int)
        d = proportions[ids].mean() - (proportions.sum() - proportions[ids].sum()) / (n - npos)
        c, t = target_n[ids].sum(), totals[ids].sum()
        r = (c / t - (target_n.sum() - c) / (totals.sum() - t)) / pooled
        p_delta += abs(d) >= abs(delta) - 1e-12
        p_roe += abs(r) >= abs(roe_delta) - 1e-12
        permutations += 1
    correction = 0 if exact else 1
    ai = rng.choice(np.flatnonzero(pos), (bootstrap, npos))
    bi = rng.choice(np.flatnonzero(neg), (bootstrap, n - npos))
    boot_delta = proportions[ai].mean(axis=1) - proportions[bi].mean(axis=1)
    ac, at = target_n[ai].sum(axis=1), totals[ai].sum(axis=1)
    bc, bt = target_n[bi].sum(axis=1), totals[bi].sum(axis=1)
    boot_pooled = (ac + bc) / (at + bt)
    with np.errstate(divide="ignore", invalid="ignore"):
        boot_roe = (ac / at) / boot_pooled
    leave_one_out = []
    for i in range(n):
        keep = np.arange(n) != i
        leave_one_out.append(proportions[keep & pos].mean() - proportions[keep & neg].mean())
    return {
        "n_positive": npos, "n_negative": n - npos,
        "mean_positive": float(a.mean()), "mean_negative": float(b.mean()),
        "median_positive": float(np.median(a)), "median_negative": float(np.median(b)),
        "delta_mean": float(delta), "delta_ci_low": float(np.quantile(boot_delta, .025)),
        "delta_ci_high": float(np.quantile(boot_delta, .975)),
        "mwu_two_sided_p": float(mannwhitneyu(a, b, alternative="two-sided", method="asymptotic").pvalue),
        "donor_permutation_p": float((p_delta + correction) / (permutations + correction)),
        "permutation_exact": exact, "donor_label_assignments": permutations,
        "roe_positive": float(roe_pos), "roe_negative": float(roe_neg),
        "roe_positive_ci_low": float(np.nanquantile(boot_roe, .025)),
        "roe_positive_ci_high": float(np.nanquantile(boot_roe, .975)),
        "roe_donor_permutation_p": float((p_roe + correction) / (permutations + correction)),
        "leave_one_donor_out_delta_min": float(min(leave_one_out)),
        "leave_one_donor_out_delta_max": float(max(leave_one_out)),
    }


def binary_cv(x, y, target, *, seeds, folds=5, epochs=200):
    torch.set_num_threads(1)
    if np.bincount(y, minlength=2).min() < folds:
        raise ValueError("Not enough donors per label for requested stratified folds")
    data = torch.tensor(x, dtype=torch.float32)
    truth = torch.tensor(y, dtype=torch.float32)
    predictions, audit = [], []
    for seed in seeds:
        oof, ablation = np.empty(len(y)), np.empty(len(y))
        for fold, (train, test) in enumerate(StratifiedKFold(folds, shuffle=True, random_state=seed).split(x, y)):
            torch.manual_seed(1000 * seed + fold)
            model = BinaryCountBagMIL(x.shape[1])
            optimizer = torch.optim.AdamW(model.parameters(), lr=.005, weight_decay=.01)
            model.train()
            for _ in range(epochs):
                optimizer.zero_grad()
                loss = nn.functional.binary_cross_entropy_with_logits(model(data[train]), truth[train])
                loss.backward()
                optimizer.step()
            model.eval()
            with torch.no_grad():
                oof[test] = model(data[test]).sigmoid().numpy()
                removed = data[test].clone()
                removed[:, target] = 0
                usable = removed.sum(dim=1) > 0
                ablation[test] = np.nan
                ablation[test[usable.numpy()]] = model(removed[usable]).sigmoid().numpy()
            audit.append({"seed": seed, "fold": fold, "train": train.tolist(), "test": test.tolist()})
        predictions.append((oof, ablation))
    probability = np.mean([p[0] for p in predictions], axis=0)
    removed_probability = np.mean([p[1] for p in predictions], axis=0)
    metrics = {"oof_auc": float(roc_auc_score(y, probability)),
               "oof_balanced_accuracy": float(balanced_accuracy_score(y, probability >= .5)),
               "seed_auc": [float(roc_auc_score(y, p[0])) for p in predictions]}
    return probability, removed_probability, metrics, audit


def _permuted_cv(task):
    index, x, y, target, cfg, directory, fingerprint = task
    path = Path(directory) / f"permutation_{index:04d}.json"
    if path.exists():
        result = json.loads(path.read_text())
        if result["fingerprint"] != fingerprint:
            raise ValueError("MIL checkpoint input/settings mismatch")
        return result
    yp = np.random.default_rng(cfg["seed"] + index + 1).permutation(y)
    _, _, metrics, _ = binary_cv(x, yp, target, seeds=cfg["mil_seeds"],
                                folds=cfg["folds"], epochs=cfg["mil_epochs"])
    result = {"permutation": index, "fingerprint": fingerprint, **metrics}
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(result) + "\n")
    tmp.replace(path)
    return result


def run_comparison(cfg, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    summary = []
    for cohort in cfg["cohorts"]:
        destination = output / cohort["name"]
        destination.mkdir(exist_ok=True)
        counts = pd.read_csv(cohort["counts"], sep="\t", index_col=0)
        metadata = pd.read_csv(cohort["metadata"], sep="\t", index_col=0)
        x, y, metadata = validate_inputs(counts, metadata)
        target = counts.columns.get_loc(cohort["target"])
        if x[:, target].sum() == 0:
            raise ValueError("Target has zero cells; Ro/e is undefined")
        stats = donor_tests(x, y, target, bootstrap=cfg["bootstrap"], seed=cfg["seed"])
        print(json.dumps({"cohort": cohort["name"], **stats}), flush=True)
        pd.DataFrame([stats]).to_csv(destination / "donor_tests.tsv", sep="\t", index=False)
        metrics = {}
        if cfg.get("run_mil", True):
            prob, removed, metrics, audit = binary_cv(x, y, target, seeds=cfg["mil_seeds"],
                                                    folds=cfg["folds"], epochs=cfg["mil_epochs"])
            predictions = metadata.copy()
            predictions["target_cells"] = x[:, target]
            predictions["dc_cells"] = x.sum(axis=1)
            predictions["target_fraction"] = x[:, target] / x.sum(axis=1)
            predictions["oof_probability"] = prob
            predictions["oof_probability_target_removed"] = removed
            predictions["oof_target_ablation_delta"] = prob - removed
            predictions.to_csv(destination / "donor_predictions.tsv", sep="\t")
            for fold in audit:
                for part in ("train", "test"):
                    fold[part] = metadata.donor_id.iloc[fold[part]].tolist()
            (destination / "fold_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
            checkpoint = destination / "mil_checkpoints"
            checkpoint.mkdir(exist_ok=True)
            signature = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode() +
                Path(cohort["counts"]).read_bytes() + Path(cohort["metadata"]).read_bytes() +
                Path(__file__).read_bytes()).hexdigest()
            tasks = [(i, x, y, target, cfg, str(checkpoint), signature) for i in range(cfg["permutations"])]
            null = []
            with ProcessPoolExecutor(max_workers=cfg["workers"], mp_context=mp.get_context("spawn")) as pool:
                futures = [pool.submit(_permuted_cv, task) for task in tasks]
                for future in as_completed(futures):
                    null.append(future.result())
                    if len(null) % 10 == 0:
                        print(f'{cohort["name"]} MIL label permutations {len(null)}/{len(tasks)}', flush=True)
            null.sort(key=lambda row: row["permutation"])
            pd.DataFrame(null).to_csv(destination / "mil_permutation_null.tsv", sep="\t", index=False)
            metrics["permutation_p_auc_greater"] = float((1 + sum(r["oof_auc"] >= metrics["oof_auc"] for r in null)) / (1 + len(null)))
            metrics["interpretation"] = "Held-out ACPA classification from input cell-state composition; no expression features. Not evidence of target expansion by itself."
            (destination / "mil_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
        summary.append({"tissue": cohort["name"], "target": cohort["target"], **stats,
                        **{k: v for k, v in metrics.items() if not isinstance(v, list)}})
    results = pd.DataFrame(summary)
    from .abundance_metrics import benjamini_hochberg
    for name in ("donor_permutation_p", "mwu_two_sided_p", "roe_donor_permutation_p", "permutation_p_auc_greater"):
        if name in results:
            results[name + "_bh"] = benjamini_hochberg(results[name].to_numpy())
    results.to_csv(output / "comparison_summary.tsv", sep="\t", index=False)
    manifest = {"method": "donor tests + Ro/e" + (" + scSurvival-style binary composition MIL" if cfg.get("run_mil", True) else ""),
                "settings": cfg, "source_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "inputs": [{"path": c[k], "sha256": hashlib.sha256(Path(c[k]).read_bytes()).hexdigest()}
                           for c in cfg["cohorts"] for k in ("counts", "metadata")],
                "versions": {"numpy": np.__version__, "pandas": pd.__version__, "torch": torch.__version__},
                "review_required": True,
                "warnings": ["Subtype composition is relative, not absolute abundance.",
                    "All p values retained; BH correction across tissues per prespecified method, no best-method selection.",
                    "MIL uses one-hot identities and binary cross-entropy; not the original expression/VAE/Cox model.",
                    "Target ablation changes relative composition; it is a prediction sensitivity, not causal attribution.",
                    "Small cohort; estimates and internal cross-validation do not establish external prediction performance."]}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return results
