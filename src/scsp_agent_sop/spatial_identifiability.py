"""Reference identifiability at measurement depth, from a v2 package's holdout cells (F12).

Pseudo-spots mix never-fitted holdout cells: each spot has ``dominant_cells`` of one component
plus ``other_cells`` drawn from the other components, binomially thinned to a log-normal depth.
They are fitted through the same ``run_rctd_api`` as production. Truth is the RNA origin of each
spot's counts, not a cell count. Reported per component: r(true, estimate), dominant recall and
top confusion; per pair: argmax confusion, mass leak and error correlation.

scRNA holdout cells carry no platform effect, so this is an optimistic bound: a pair flagged
here is not separable at this depth; an unflagged pair may still not be.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from scoop.io import atomic_json
from .rctd_adapter import run_rctd_api
from .storage import write_table

LABEL = "identifiability"


def build_pseudospots(holdout, genes, *, spots_per_component=40, dominant_cells=3, other_cells=3,
                      depth_median=4000.0, depth_sdlog=0.6, seed=20261031):
    """Return (count view on ``genes``, truth RNA fractions, dominant component per spot)."""
    rng = np.random.default_rng(seed)
    x = sparse.csr_matrix(holdout.X).astype(np.int64)
    component = holdout.obs.component.astype(str).to_numpy()
    usable = sorted(set(component))
    if len(usable) < 2:
        raise ValueError("Pseudo-spots need holdout cells from at least two components")
    pools = {t: np.flatnonzero(component == t) for t in usable}
    index = {t: i for i, t in enumerate(usable)}
    n_spots = spots_per_component * len(usable)
    member, counts, dominant = np.zeros((n_spots, len(usable))), [], np.empty(n_spots, int)
    for s in range(n_spots):
        dom = usable[s % len(usable)]
        others = rng.choice([t for t in usable if t != dom], other_cells, replace=True)
        cells = [(dom, c) for c in rng.choice(pools[dom], dominant_cells, replace=True)]
        cells += [(t, rng.choice(pools[t])) for t in others]
        rows = x[[c for _, c in cells]]
        per_cell = np.asarray(rows.sum(axis=1)).ravel().astype(float)
        keep = min(1.0, rng.lognormal(np.log(depth_median), depth_sdlog) / per_cell.sum())
        counts.append(sparse.csr_matrix(rng.binomial(np.asarray(rows.sum(axis=0)).ravel(), keep)))
        for (t, _), mass in zip(cells, per_cell):
            member[s, index[t]] += mass
        dominant[s] = index[dom]
    y = sparse.vstack(counts).tocsr()
    names = [f"ps{i}" for i in range(n_spots)]
    import anndata as ad

    columns = pd.Index(holdout.var_names.astype(str)).get_indexer(genes)
    view = ad.AnnData(y[:, columns].astype(np.float64),
                      obs=pd.DataFrame({"native_full_exposure": np.asarray(y.sum(axis=1)).ravel()}, index=names),
                      var=pd.DataFrame({"is_gene_expression": True}, index=pd.Index(genes)))
    truth = pd.DataFrame(member / member.sum(axis=1, keepdims=True), index=names, columns=usable)
    return view, truth, np.asarray(usable)[dominant]


def score(estimate: pd.DataFrame, truth: pd.DataFrame, dominant, *, confusion_flag=0.2, leak_flag=0.15, errcorr_flag=-0.4):
    """Per-component recovery and per-pair separability on fitted pseudo-spots."""
    types = list(truth.columns)
    E, T = estimate[types].to_numpy(), truth.to_numpy()
    D = pd.Index(types).get_indexer(dominant)
    est_dom = E.argmax(axis=1)
    K = len(types)

    def r(a, b):
        return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else np.nan

    conf, leak = np.zeros((K, K)), np.zeros((K, K))
    per = []
    for k, name in enumerate(types):
        m = D == k
        conf[k] = np.bincount(est_dom[m], minlength=K) / max(m.sum(), 1)
        leak[k] = (E[m] - T[m]).mean(axis=0) if m.any() else 0
        per.append({"component": name, "n_dominant_spots": int(m.sum()), "r_true_vs_est": r(T[:, k], E[:, k]),
                    "dominant_recall": float(conf[k, k]) if m.any() else np.nan,
                    "mean_est_when_absent": float(E[T[:, k] == 0, k].mean()) if (T[:, k] == 0).any() else np.nan,
                    "top_confusion": types[int(np.argmax(other))] if (other := np.where(np.arange(K) == k, 0, conf[k])).max() > 0 else ""})
    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.corrcoef((E - T).T)
    pairs, partners = [], {t: [] for t in types}
    for i in range(K):
        for j in range(i + 1, K):
            mean_conf, max_leak = (conf[i, j] + conf[j, i]) / 2, max(leak[i, j], leak[j, i])
            flag = bool(mean_conf >= confusion_flag or max_leak >= leak_flag or C[i, j] <= errcorr_flag)
            pairs.append({"A": types[i], "B": types[j], "confusion_A_to_B": conf[i, j], "confusion_B_to_A": conf[j, i],
                          "mean_confusion": mean_conf, "leak_A_to_B": leak[i, j], "leak_B_to_A": leak[j, i],
                          "error_corr": C[i, j], "r_merged_true_vs_est": r(T[:, i] + T[:, j], E[:, i] + E[:, j]),
                          "non_separable": flag})
            # Confusable partners: flagged, or a weaker leak / error coupling (project rule).
            if flag or max_leak >= 0.05 or C[i, j] <= -0.2:
                partners[types[i]].append(types[j])
                partners[types[j]].append(types[i])
    per = pd.DataFrame(per)
    per["confusable_partners"] = per.component.map(lambda c: ";".join(sorted(partners[c])))
    per["identifiability_verdict"] = np.select([per.r_true_vs_est >= 0.7, per.r_true_vs_est >= 0.4],
                                               ["identifiable", "weakly_identifiable"], "not_identifiable_fine")
    pairs = pd.DataFrame(pairs).sort_values(["non_separable", "mean_confusion"], ascending=False)
    return per, pairs


def run_identifiability(reference_dir, output_dir, *, rctd_options=None, run_purpose="pilot", preset="standard",
                        spots_per_component=40, dominant_cells=3, other_cells=3, depth_median=4000.0, depth_sdlog=0.6,
                        confusion_flag=0.2, leak_flag=0.15, errcorr_flag=-0.4, seed=20261031) -> dict:
    """Fit holdout pseudo-spots of a v2 package and write component/pair identifiability tables.

    ``component_identifiability.parquet`` is the table ``spatial.rctd.identifiability_table`` joins
    into each fit's component summary.
    """
    reference, out = Path(reference_dir), Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = reference / "reference_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != "scoop.reference-profiles.v2":
        raise ValueError("The identifiability harness needs a scoop.reference-profiles.v2 package with holdout cells")
    import anndata as ad

    profile_path = reference / manifest["profiles_path"]
    with np.load(profile_path, allow_pickle=False) as z:
        genes, types = z["gene_names"].astype(str).tolist(), z["cell_type_names"].astype(str).tolist()
    holdout = ad.read_h5ad(reference / manifest["holdout_rule"]["path"])
    view, truth, dominant = build_pseudospots(holdout, genes, spots_per_component=spots_per_component,
        dominant_cells=dominant_cells, other_cells=other_cells, depth_median=depth_median, depth_sdlog=depth_sdlog, seed=seed)
    view_path = out / "pseudospots_view.h5ad"
    view.write_h5ad(view_path)
    write_table(truth.rename_axis("spatial_unit_id").reset_index(), out / "pseudospots_truth_rna_fraction.parquet")
    heldout = reference / manifest["heldout_validation_genes"]["path"]
    exclude = pd.read_csv(heldout, sep="\t").gene.tolist() if heldout.is_file() else []
    options = {"reference_profiles": str(profile_path), "reference_manifest": str(manifest_path),
               "spatial_counts_layer": "X", "spatial_exposure_key": "native_full_exposure", "umi_min": 100,
               "preset": preset, "device": "cpu", "seed": seed, "exclude_genes": exclude,
               **(rctd_options or {}), "run_purpose": run_purpose, "output_label": LABEL}
    cfg = {"qc": {"counts_layer": "X", "gene_expression_key": "is_gene_expression"}, "spatial": {"rctd": options}}
    info = run_rctd_api(view_path, None, out / "rctd_fit" / f"rctd_{LABEL}.h5ad", cfg, "full")
    weights = pd.read_parquet(out / "rctd_fit" / f"{LABEL}_weights_visualization.parquet").set_index("spatial_unit_id")
    fitted = weights.notna().all(axis=1).to_numpy()
    estimate = weights.loc[fitted].reindex(columns=types)
    per, pairs = score(estimate, truth.loc[fitted].reindex(columns=types, fill_value=0.0), dominant[fitted],
                       confusion_flag=confusion_flag, leak_flag=leak_flag, errcorr_flag=errcorr_flag)
    declared = pd.DataFrame(manifest["components"]).rename(columns={"name": "component"})
    per = declared[["component", "review_status", "holdout_cells", "holdout_donors", "holdout_status"]].merge(per, on="component")
    # A component without holdout cells never dominates a spot: its recovery is not assessed.
    per.loc[per.n_dominant_spots == 0, "identifiability_verdict"] = "not_assessable"
    write_table(per, out / "component_identifiability.parquet")
    write_table(pairs, out / "pair_identifiability.parquet")
    summary = {"reference_manifest": str(manifest_path.resolve()), "profiles_sha256": manifest["profiles_sha256"],
               "n_components": len(types), "n_pseudospots": len(view), "n_fitted": int(fitted.sum()),
               "design": {"dominant_cells": dominant_cells, "other_cells": other_cells, "depth_median": depth_median,
                          "depth_sdlog": depth_sdlog, "spots_per_component": spots_per_component, "preset": preset,
                          "heldout_genes_excluded": len(exclude), "seed": seed},
               "flag_rule": {"mean_confusion_ge": confusion_flag, "leak_ge": leak_flag, "error_corr_le": errcorr_flag},
               "median_r_true_vs_est": float(per.r_true_vs_est.median()),
               "median_dominant_recall": float(per.dominant_recall.median()),
               "verdicts": per.identifiability_verdict.value_counts().to_dict(),
               "n_non_separable_pairs": int(pairs.non_separable.sum()),
               "heldout_gene_validation": info["heldout_gene_validation"], "sigma": info["sigma"],
               "fit_identity": info["fit_identity"],
               "caveat": "scRNA holdout pseudo-spots without a platform effect: optimistic bound; truth is RNA origin"}
    atomic_json(out / "identifiability_summary.json", summary)
    return summary
