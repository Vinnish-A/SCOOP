from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import optimize, special, stats


@dataclass(frozen=True)
class FastDEResult:
    table: pd.DataFrame
    design: pd.DataFrame
    size_factors: pd.Series
    dispersions: pd.Series
    manifest: dict


def run_deseq2_wald(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    condition_col: str,
    ctrl_group: str,
    test_group: str,
    min_total_count: int = 10,
    min_samples_per_group: int = 2,
    min_dispersion: float = 1e-8,
    max_iter: int = 50,
    beta_tol: float = 1e-8,
) -> FastDEResult:
    """DESeq2-like pseudobulk differential expression with NB Wald tests.

    The implementation mirrors the main DESeq2 concepts used for a two-group
    pseudobulk contrast: median-ratio size factors, gene-wise negative-binomial
    dispersion estimates, parametric dispersion trend/MAP shrinkage, log-link
    GLM with size-factor offset, and a Wald test for the condition coefficient.
    It intentionally does not perform LFC shrinkage, matching DESeq2's default
    DESeq(..., betaPrior=FALSE) + results(...) path.
    """

    start = perf_counter()
    if condition_col not in metadata.columns:
        raise KeyError(f"condition_col {condition_col!r} is not present in metadata")
    meta = metadata[metadata[condition_col].astype(str).isin([str(ctrl_group), str(test_group)])].copy()
    counts = counts.loc[meta.index]
    condition = meta[condition_col].astype(str)
    ctrl_mask = condition == str(ctrl_group)
    test_mask = condition == str(test_group)
    if int(ctrl_mask.sum()) < min_samples_per_group or int(test_mask.sum()) < min_samples_per_group:
        raise ValueError("not enough pseudobulk samples in one or both groups")

    y_all = counts.to_numpy(dtype=np.float64, copy=False)
    keep = y_all.sum(axis=0) >= min_total_count
    y = y_all[:, keep]
    genes = counts.columns[keep].astype(str)
    size_factors = estimate_size_factors(y)
    norm = y / size_factors[:, None]
    base_mean = norm.mean(axis=0)
    x = test_mask.astype(float).to_numpy()
    design_np = np.column_stack([np.ones_like(x), x])
    initial_disp = estimate_initial_dispersions(norm, design_np, min_dispersion=min_dispersion)
    beta_init, _ = fit_nb_glm_binary_condition(
        y,
        x=x,
        size_factors=size_factors,
        dispersions=initial_disp,
        max_iter=max_iter,
        beta_tol=beta_tol,
    )
    mu_init = fitted_mu_binary(beta_init, x=x, size_factors=size_factors, min_mu=0.5)
    disp_gene = estimate_gene_wise_dispersions_mle(
        y,
        mu=mu_init,
        x=x,
        min_dispersion=min_dispersion,
        max_dispersion=max(10.0, float(y.shape[0])),
    )
    disp_fit, trend_info = fit_parametric_dispersion_trend(base_mean, disp_gene, min_dispersion=min_dispersion)
    disp_prior_var = estimate_dispersion_prior_var(
        disp_gene,
        disp_fit,
        n_samples=y.shape[0],
        n_coefficients=design_np.shape[1],
        min_dispersion=min_dispersion,
        fallback=0.25,
    )
    disp_map = estimate_map_dispersions(
        y,
        mu=mu_init,
        x=x,
        disp_gene=disp_gene,
        disp_fit=disp_fit,
        prior_var=disp_prior_var,
        min_dispersion=min_dispersion,
        max_dispersion=max(10.0, float(y.shape[0])),
    )
    var_log_disp_ests = float(np.nanvar(np.log(np.maximum(disp_gene, min_dispersion)) - np.log(np.maximum(disp_fit, min_dispersion)), ddof=1))
    disp_outlier = np.log(np.maximum(disp_gene, min_dispersion)) > (
        np.log(np.maximum(disp_fit, min_dispersion)) + 2.0 * np.sqrt(max(var_log_disp_ests, 0.0))
    )
    dispersions = np.where(disp_outlier, disp_gene, disp_map)
    beta, se = fit_nb_glm_binary_condition(
        y,
        x=x,
        size_factors=size_factors,
        dispersions=dispersions,
        max_iter=max_iter,
        beta_tol=beta_tol,
    )
    beta_condition = beta[1]
    se_condition = se
    stat = np.divide(beta_condition, se_condition, out=np.zeros_like(beta_condition), where=se_condition > 0)
    pvalue = 2.0 * stats.norm.sf(np.abs(stat))
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    padj = benjamini_hochberg(pvalue)
    log2fc = beta_condition / np.log(2.0)
    lfc_se = se_condition / np.log(2.0)

    table = pd.DataFrame(
        {
            "gene": genes,
            "baseMean": base_mean,
            "log2FoldChange": log2fc,
            "lfcSE": lfc_se,
            "stat": stat,
            "pvalue": pvalue,
            "padj": padj,
            "dispersion": dispersions,
            "dispGeneEst": disp_gene,
            "dispFit": disp_fit,
            "dispMAP": disp_map,
            "dispOutlier": disp_outlier,
        }
    ).sort_values(["pvalue", "gene"], kind="mergesort").reset_index(drop=True)
    design = pd.DataFrame(
        {
            "Intercept": 1.0,
            f"{condition_col}_{test_group}_vs_{ctrl_group}": x,
        },
        index=meta.index,
    )
    sf = pd.Series(size_factors, index=meta.index, name="size_factor")
    disp = pd.Series(dispersions, index=genes, name="dispersion")
    manifest = {
        "schema_version": "fastde.deseq2_wald.v1",
        "method": "fastde_deseq2_wald",
        "condition_col": condition_col,
        "ctrl_group": str(ctrl_group),
        "test_group": str(test_group),
        "n_samples": int(len(meta)),
        "n_ctrl": int(ctrl_mask.sum()),
        "n_test": int(test_mask.sum()),
        "n_input_genes": int(counts.shape[1]),
        "n_tested_genes": int(table.shape[0]),
        "min_total_count": int(min_total_count),
        "max_iter": int(max_iter),
        "beta_tol": float(beta_tol),
        "fit_type": "parametric",
        "dispersion_trend": trend_info,
        "dispersion_prior_var": float(disp_prior_var),
        "seconds": round(perf_counter() - start, 6),
        "limitations": "DESeq2-like two-group Wald path; no LFC shrinkage or complex design support yet",
    }
    return FastDEResult(table=table, design=design, size_factors=sf, dispersions=disp, manifest=manifest)


def estimate_size_factors(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.float64)
    positive_all = np.all(counts > 0, axis=0)
    if np.any(positive_all):
        logs = np.log(counts[:, positive_all])
        geo = np.exp(logs.mean(axis=0))
    else:
        positive = counts > 0
        logs = np.where(positive, np.log(np.maximum(counts, 1.0)), np.nan)
        geo = np.exp(np.nanmean(logs, axis=0))
        positive_all = np.isfinite(geo) & (geo > 0)
        logs = np.log(np.maximum(counts[:, positive_all], 1.0))
        geo = geo[positive_all]
    ratios = counts[:, positive_all] / geo[None, :]
    ratios = np.where(ratios > 0, ratios, np.nan)
    size_factors = np.nanmedian(ratios, axis=1)
    if not np.all(np.isfinite(size_factors)) or np.any(size_factors <= 0):
        lib = counts.sum(axis=1)
        size_factors = lib / np.exp(np.mean(np.log(np.maximum(lib, 1.0))))
    size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))
    return size_factors


def estimate_gene_wise_dispersions(norm_counts: np.ndarray, *, min_dispersion: float) -> np.ndarray:
    mean = norm_counts.mean(axis=0)
    var = norm_counts.var(axis=0, ddof=1)
    disp = (var - mean) / np.maximum(mean**2, np.finfo(float).tiny)
    disp = np.where(np.isfinite(disp), disp, min_dispersion)
    return np.maximum(disp, min_dispersion)


def estimate_initial_dispersions(norm_counts: np.ndarray, design: np.ndarray, *, min_dispersion: float) -> np.ndarray:
    """DESeq2-style starting values: min(rough design-aware, moments)."""

    y = np.asarray(norm_counts, dtype=np.float64)
    x = np.asarray(design, dtype=np.float64)
    beta = np.linalg.pinv(x) @ y
    mu = np.maximum(x @ beta, 1.0)
    m, p = x.shape
    denom = max(m - p, 1)
    rough = (((y - mu) ** 2 - mu) / np.maximum(mu**2, np.finfo(float).tiny)).sum(axis=0) / denom
    moments = estimate_gene_wise_dispersions(y, min_dispersion=min_dispersion)
    disp = np.minimum(np.maximum(rough, 0.0), moments)
    disp = np.where(np.isfinite(disp), disp, min_dispersion)
    return np.maximum(disp, min_dispersion)


def estimate_gene_wise_dispersions_mle(
    counts: np.ndarray,
    *,
    mu: np.ndarray,
    x: np.ndarray,
    min_dispersion: float,
    max_dispersion: float,
) -> np.ndarray:
    y = np.asarray(counts, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    lower = np.log(min_dispersion / 10.0)
    upper = np.log(max_dispersion)
    out = np.empty(y.shape[1], dtype=np.float64)
    for j in range(y.shape[1]):
        res = optimize.minimize_scalar(
            lambda log_alpha: -_cox_reid_loglik_alpha(y[:, j], mu[:, j], x, np.exp(log_alpha)),
            bounds=(lower, upper),
            method="bounded",
            options={"xatol": 1e-4},
        )
        out[j] = np.exp(res.x) if res.success else min_dispersion
    return np.clip(out, min_dispersion, max_dispersion)


def fit_parametric_dispersion_trend(base_mean: np.ndarray, disp_gene: np.ndarray, *, min_dispersion: float) -> tuple[np.ndarray, dict]:
    means = np.asarray(base_mean, dtype=np.float64)
    disps = np.asarray(disp_gene, dtype=np.float64)
    use = np.isfinite(means) & np.isfinite(disps) & (means > 0) & (disps > 100.0 * min_dispersion)
    if int(use.sum()) < 10:
        mean_disp = float(np.nanmean(disps[disps > 10.0 * min_dispersion]))
        fit = np.full_like(disps, max(mean_disp, min_dispersion))
        return fit, {"type": "mean", "mean": float(fit[0])}

    coefs = np.array([0.1, 1.0], dtype=np.float64)
    ok = use.copy()
    converged = False
    for _ in range(20):
        m = means[ok]
        d = disps[ok]
        if len(m) < 10:
            break

        def residual(log_coefs: np.ndarray) -> np.ndarray:
            a, b = np.exp(log_coefs)
            pred = a + b / m
            # Gamma(identity)-like fit: relative residuals approximate the
            # robust filtering used by DESeq2's parametricDispersionFit.
            return (d - pred) / np.maximum(pred, min_dispersion)

        trial = optimize.least_squares(residual, np.log(coefs), max_nfev=2000)
        if not trial.success:
            break
        new_coefs = np.exp(trial.x)
        pred_all = new_coefs[0] + new_coefs[1] / means
        rel = disps / np.maximum(pred_all, min_dispersion)
        new_ok = use & (rel > 1e-4) & (rel < 15.0)
        converged = np.sum(np.log(new_coefs / coefs) ** 2) < 1e-6
        coefs = new_coefs
        ok = new_ok
        if converged:
            break

    if not np.all(np.isfinite(coefs)) or np.any(coefs <= 0):
        mean_disp = float(np.nanmean(disps[use]))
        fit = np.full_like(disps, max(mean_disp, min_dispersion))
        return fit, {"type": "mean_fallback", "mean": float(fit[0])}
    fit = np.maximum(coefs[0] + coefs[1] / np.maximum(means, np.finfo(float).tiny), min_dispersion)
    return fit, {"type": "parametric", "asymptDisp": float(coefs[0]), "extraPois": float(coefs[1]), "converged": bool(converged)}


def estimate_dispersion_prior_var(
    disp_gene: np.ndarray,
    disp_fit: np.ndarray,
    *,
    n_samples: int,
    n_coefficients: int,
    min_dispersion: float,
    fallback: float,
) -> float:
    above = disp_gene >= min_dispersion * 100.0
    resid = np.log(np.maximum(disp_gene[above], min_dispersion)) - np.log(np.maximum(disp_fit[above], min_dispersion))
    resid = resid[np.isfinite(resid)]
    if resid.size < 2:
        return float(fallback)
    df = n_samples - n_coefficients
    var_log_disp_ests = float(np.var(resid, ddof=1))
    if df > 3:
        return float(max(var_log_disp_ests - special.polygamma(1, df / 2.0), fallback))
    return float(fallback)


def estimate_map_dispersions(
    counts: np.ndarray,
    *,
    mu: np.ndarray,
    x: np.ndarray,
    disp_gene: np.ndarray,
    disp_fit: np.ndarray,
    prior_var: float,
    min_dispersion: float,
    max_dispersion: float,
) -> np.ndarray:
    y = np.asarray(counts, dtype=np.float64)
    lower = np.log(min_dispersion / 10.0)
    upper = np.log(max_dispersion)
    out = np.empty(y.shape[1], dtype=np.float64)
    prior_var = max(float(prior_var), 1e-8)
    for j in range(y.shape[1]):
        init = disp_gene[j] if disp_gene[j] > 0.1 * disp_fit[j] else disp_fit[j]
        init = float(np.clip(init, min_dispersion, max_dispersion))
        prior_mean = np.log(max(disp_fit[j], min_dispersion))

        def objective(log_alpha: float) -> float:
            alpha = np.exp(log_alpha)
            prior = 0.5 * ((log_alpha - prior_mean) ** 2) / prior_var
            return -_cox_reid_loglik_alpha(y[:, j], mu[:, j], x, alpha) + prior

        lo = max(lower, np.log(init) - 4.0)
        hi = min(upper, np.log(init) + 4.0)
        res = optimize.minimize_scalar(objective, bounds=(lo, hi), method="bounded", options={"xatol": 1e-4})
        out[j] = np.exp(res.x) if res.success else init
    return np.clip(out, min_dispersion, max_dispersion)


def _cox_reid_loglik_alpha(y: np.ndarray, mu: np.ndarray, x: np.ndarray, alpha: float) -> float:
    mu = np.maximum(np.asarray(mu, dtype=np.float64), 0.5)
    y = np.asarray(y, dtype=np.float64)
    alpha = float(max(alpha, np.finfo(float).tiny))
    size = 1.0 / alpha
    ll = (
        special.gammaln(y + size)
        - special.gammaln(size)
        - special.gammaln(y + 1.0)
        + size * (np.log(size) - np.log(size + mu))
        + y * (np.log(mu) - np.log(size + mu))
    ).sum()
    weights = mu / (1.0 + alpha * mu)
    a = weights.sum()
    b = (weights * x).sum()
    c = (weights * x * x).sum()
    det = max(a * c - b * b, np.finfo(float).tiny)
    return float(ll - 0.5 * np.log(det))


def fit_nb_glm_binary_condition(
    counts: np.ndarray,
    *,
    x: np.ndarray,
    size_factors: np.ndarray,
    dispersions: np.ndarray,
    max_iter: int,
    beta_tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(counts, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    offset = np.log(size_factors)
    norm = y / size_factors[:, None]
    ctrl = x == 0
    test = x == 1
    beta = np.zeros((2, y.shape[1]), dtype=np.float64)
    beta[0] = np.log(np.maximum(norm[ctrl].mean(axis=0), 1e-8))
    beta[1] = np.log(np.maximum(norm[test].mean(axis=0), 1e-8)) - beta[0]
    alpha = dispersions[None, :]

    for _ in range(max_iter):
        eta_no_offset = beta[0][None, :] + x[:, None] * beta[1][None, :]
        eta = offset[:, None] + eta_no_offset
        eta = np.clip(eta, -30.0, 30.0)
        mu = np.maximum(np.exp(eta), 0.5)
        weights = mu / (1.0 + alpha * mu)
        z = eta_no_offset + (y - mu) / np.maximum(mu, np.finfo(float).tiny)
        a = weights.sum(axis=0)
        b = (weights * x[:, None]).sum(axis=0)
        c = (weights * (x[:, None] ** 2)).sum(axis=0)
        d = (weights * z).sum(axis=0)
        e = (weights * x[:, None] * z).sum(axis=0)
        det = a * c - b * b
        valid = np.abs(det) > np.finfo(float).tiny
        new_beta = beta.copy()
        new_beta[0, valid] = (d[valid] * c[valid] - b[valid] * e[valid]) / det[valid]
        new_beta[1, valid] = (a[valid] * e[valid] - b[valid] * d[valid]) / det[valid]
        if np.nanmax(np.abs(new_beta - beta)) < beta_tol:
            beta = new_beta
            break
        beta = new_beta

    eta_no_offset = beta[0][None, :] + x[:, None] * beta[1][None, :]
    mu = np.maximum(np.exp(np.clip(offset[:, None] + eta_no_offset, -30.0, 30.0)), 0.5)
    weights = mu / (1.0 + alpha * mu)
    a = weights.sum(axis=0)
    b = (weights * x[:, None]).sum(axis=0)
    c = (weights * (x[:, None] ** 2)).sum(axis=0)
    det = a * c - b * b
    se = np.sqrt(np.divide(a, det, out=np.full_like(a, np.inf), where=det > 0))
    return beta, se


def fitted_mu_binary(beta: np.ndarray, *, x: np.ndarray, size_factors: np.ndarray, min_mu: float) -> np.ndarray:
    offset = np.log(size_factors)
    eta = offset[:, None] + beta[0][None, :] + x[:, None] * beta[1][None, :]
    return np.maximum(np.exp(np.clip(eta, -30.0, 30.0)), min_mu)


def benjamini_hochberg(pvalue: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalue, dtype=float)
    n = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out
