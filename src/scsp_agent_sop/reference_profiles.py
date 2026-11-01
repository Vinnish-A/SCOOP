"""Measurement-aware RCTD reference package producer (``scoop.reference-profiles.v2``).

One streamed pass over native counts through ``H5ADReader``. For component k and gene g

    profile[g, k] = sum_{c in k, sample(c) measured g} x_cg / nUMI_c  /  n_measured[g, k]

where nUMI_c is the cell's full native exposure and n_measured[g, k] counts the training
cells whose sample measured g (``FeatureMeasurement``; without one, every cell measured every
gene). This is ``rctd.Reference``'s estimator, mean(x / nUMI), with unmeasured entries left out
instead of counted as zeros, over every admitted training cell (no per-type cap or subsample).

Sums are kept per (component, group), so merged packages re-aggregate exactly from
``component_sums.npz`` without another pass. A per-component cell holdout, chosen by a declared
rule, is excluded from the sums and saved for identifiability tests.

``native_exposure`` and ``component_review_status`` are the exposure and review rules the RCTD
adapter applies to its own H5AD inputs; the package records the basis of each.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse

from scoop.io import atomic_json, sha256_file
from . import __version__
from .assays import gene_expression_mask
from .h5ad_io import H5ADReader
from .measurement import var_sha256

SCHEMA = "scoop.reference-profiles.v2"
UNIT = "mean_counts_per_full_native_exposure"
# In increasing order of caution; a component takes the most cautious status among its cells.
REVIEW_STATUSES = ("reviewed", "low_confidence", "unresolved", "evidence_pending")
# SCOOP's committed broad_/fine_annotation_status (annotation_decision.committer) read as a review
# status; a level never decided (unassigned) or not applicable is not reviewed.
COMMITTED_STATUS = {"accepted": "reviewed", "unresolved": "unresolved", "review_required": "evidence_pending",
                    "unassigned": "evidence_pending", "not_applicable": "evidence_pending"}
ANNOTATION_LEVELS = ("broad", "fine")
RCTD_FC, RCTD_EXPR = 0.75, 2e-4  # rctd.Reference.get_de_genes defaults
TRUNCATED = "gene_axis_row_sum_with_unmeasured_genes"


def native_exposure(obs, total, exposure_key=None):
    """Full native exposure per row and its basis.

    A declared ``exposure_key`` must be integer depth >= the row's gene counts ``total`` (None leaves
    that comparison to a caller that streams the counts). Otherwise the gene-axis row sum is the
    exposure, and it must equal a bound ``native_full_gene_UMI`` wherever that is declared.
    """
    if exposure_key:
        exposure = pd.to_numeric(obs[exposure_key], errors="raise").to_numpy(float)
        if (not np.isfinite(exposure).all() or not np.equal(exposure, np.rint(exposure)).all()
                or (total is not None and (exposure < total).any())):
            raise ValueError("Declared full native exposure must be finite integer depth >= measured count mass")
        return exposure, f"declared:{exposure_key}"
    if "native_full_gene_UMI" in obs:
        declared = pd.to_numeric(obs.native_full_gene_UMI, errors="raise").to_numpy(float)
        measured = np.isfinite(declared)
        if not np.array_equal(declared[measured], total[measured]):
            raise ValueError("Native count depth differs from the bound observation registry")
        if measured.all():
            return total, "gene_axis_row_sum_matches_native_full_gene_UMI"
    return total, "gene_axis_row_sum"


def component_review_status(obs, k_cell, n_components, *, level="fine", review_status_key=None):
    """Each component's review status, the most cautious over its cells and over every review source.

    Sources: SCOOP's committed ``{level}_annotation_status`` at the label's annotation level, and an
    explicit per-cell ``review_status_key`` (review statuses or committed values). Without either, a
    component is ``evidence_pending``. Returns the statuses and the sources read.
    """
    if level not in ANNOTATION_LEVELS:
        raise ValueError(f"annotation_level must be one of {ANNOTATION_LEVELS}")
    committed = f"{level}_annotation_status"
    sources = [committed] if committed in obs else []
    if review_status_key and review_status_key != committed:
        sources.append(review_status_key)
    member = np.asarray(k_cell) >= 0
    rank = np.full(int(member.sum()), 0 if sources else REVIEW_STATUSES.index("evidence_pending"))
    for key in sources:
        cell = pd.Series(obs[key].to_numpy(object)[member]).map(
            lambda v: "evidence_pending" if pd.isna(v) else COMMITTED_STATUS.get(str(v), str(v)))
        unknown = set(cell) - set(REVIEW_STATUSES)
        if unknown:
            raise ValueError(f"Unknown review status in obs[{key!r}]: {sorted(unknown)}")
        rank = np.maximum(rank, cell.map(REVIEW_STATUSES.index).to_numpy())
    worst = pd.Series(rank).groupby(np.asarray(k_cell)[member]).max()
    return [REVIEW_STATUSES[worst.get(k, len(REVIEW_STATUSES) - 1)] for k in range(n_components)], sources


def registered_measurement(path) -> Path | None:
    """The ``feature_measurement`` artifact registered in an H5AD's uns, read without loading it."""
    from anndata.io import read_elem

    with H5ADReader(path) as reader:
        node = reader.file.get("uns/file_registry/artifacts/feature_measurement/path")
        return Path(read_elem(node)) if node is not None else None


def _counts(reader, rows, genes, layer):
    x = reader.read_rows(rows, genes, layer=layer)
    x = sparse.csr_matrix(x)
    if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.equal(x.data, np.rint(x.data)).all():
        raise ValueError("Reference counts must be finite nonnegative integers, not normalized expression")
    return x


def _row_sums(path, layer, genes, rows, chunk):
    with H5ADReader(path) as reader:
        return np.concatenate([np.asarray(_counts(reader, rows[a:a + chunk], genes, layer).sum(axis=1)).ravel()
                               for a in range(0, len(rows), chunk)] or [np.empty(0)])


def _stream(path, layer, genes, rows, combo, exposure, n_combo, holdout, chunk):
    """Sum x/nUMI per (component, group) over ``rows``; collect holdout rows' counts."""
    sums = np.zeros((n_combo, len(genes)))
    held = []
    with H5ADReader(path) as reader:
        for a in range(0, len(rows), chunk):
            part = rows[a:a + chunk]
            x = _counts(reader, part, genes, layer)
            if (exposure[part] < np.asarray(x.sum(axis=1)).ravel()).any():
                raise ValueError("Declared native exposure is below a cell's measured gene counts")
            hold = np.flatnonzero(holdout[part])
            if len(hold):
                held.append((part[hold], x[hold]))
            fit = np.flatnonzero(combo[part] >= 0)
            assign = sparse.csr_matrix((1.0 / exposure[part[fit]], (combo[part[fit]], np.arange(len(fit)))),
                                       shape=(n_combo, len(fit)))
            product = (assign @ x[fit].astype(np.float64)).tocsr()
            product.sum_duplicates()
            sums[np.repeat(np.arange(n_combo), np.diff(product.indptr)), product.indices] += product.data
    return sums, held


def _map_jobs(function, jobs, workers):
    if workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(workers) as pool:
            return list(pool.map(function, *zip(*jobs)))
    return [function(*job) for job in jobs]


def _per_component(values, k_cell, names, field):
    """One value per component from a per-cell column; several values are an error."""
    found = pd.Series(values).groupby(k_cell).agg(lambda s: sorted(set(s.dropna().astype(str))))
    mixed = [names[k] for k, v in found.items() if len(v) > 1]
    if mixed:
        raise ValueError(f"Components map to several {field} values: {mixed[:5]}")
    return [found[k][0] if found[k] else None for k in range(len(names))]


def build_reference(counts_h5ad, output_dir, *, label_key: str, donor_key: str, labels: pd.DataFrame | None = None,
                    layer: str | None = None, status_key: str = "annotation_status",
                    admit_status=("accepted", "labeled"), reviewed_key: str = "annotation_reviewed",
                    technical_key: str | None = None, require_true=(), review_status_key: str | None = None,
                    annotation_level: str = "fine",
                    parent_key: str | None = None, compartment_key: str | None = None,
                    components: pd.DataFrame | None = None, group_key: str | None = None,
                    chemistry_key: str | None = None, exposure_key: str | None = None, measurement=None,
                    min_umi: int = 100, min_cells: int = 25, min_donors: int = 2,
                    min_measured_fraction: float = 0.5, min_measured_cells: int = 25,
                    holdout_fraction: float = 0.1, holdout_max_cells: int = 200, holdout_min_cells: int = 20,
                    heldout_genes_per_component: int = 3, heldout_random_genes: int = 50,
                    workers: int = 1, chunk_rows: int = 50000, seed: int = 20261030, cfg=None) -> dict:
    """Write a v2 reference package to ``output_dir`` and return its manifest.

    Cells enter only with ``status_key`` in ``admit_status``, ``reviewed_key`` true (when present),
    every ``require_true`` column true and exposure >= ``min_umi``. ``excluded_technical`` cells (and
    ``technical_key`` cells) are never admitted. Review status comes from ``component_review_status``
    at ``annotation_level``. ``components`` (columns component, and any of parent, compartment,
    review_status) overrides per-cell parent/compartment; its review_status can only add caution.
    """
    started = perf_counter()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    path = Path(counts_h5ad)
    with H5ADReader(path) as reader:
        obs, var_names = reader.obs, reader.var_names.astype(str)
        genes = np.flatnonzero(gene_expression_mask(reader, cfg or {"run": {"assay": "scrna"}}))
    if not obs.index.is_unique or not var_names.is_unique:
        raise ValueError("Reference observation/gene IDs must be unique; reconcile native namespaces/duplicate genes explicitly")
    if labels is not None:
        missing = ~obs.index.isin(labels.index)
        if missing.any():
            raise ValueError(f"Labels table does not cover {int(missing.sum())} observations")
        obs = obs.drop(columns=[c for c in labels.columns if c in obs]).join(labels)
    if measurement is not None:
        measurement.check(SimpleNamespace(obs=obs, var_names=var_names))
    n = len(obs)

    # ---- admission
    status = obs[status_key].astype(str)
    technical = status.eq("excluded_technical").to_numpy()
    if technical_key:
        if not obs[technical_key].dropna().isin([True, False, 0, 1]).all():
            raise ValueError(f"obs[{technical_key!r}] must be boolean")
        technical |= obs[technical_key].eq(True).fillna(False).to_numpy(bool)
    admitted = status.isin(admit_status).to_numpy() & ~technical
    for key in ([reviewed_key] if reviewed_key in obs else []) + list(require_true):
        admitted &= obs[key].eq(True).fillna(False).to_numpy(bool)
    label = obs[label_key]
    if label[admitted].isna().any() or label[admitted].astype(str).str.strip().eq("").any():
        raise ValueError("Admitted reference cells need nonmissing component labels")
    rows = np.flatnonzero(admitted)
    total = None  # a declared exposure is checked against each row's counts in the streamed pass
    if not exposure_key:
        jobs = [(path, layer, genes, block, chunk_rows) for block in np.array_split(rows, max(1, 2 * workers)) if len(block)]
        total = np.concatenate(_map_jobs(_row_sums, jobs, workers) or [np.empty(0)])
    exposure = np.zeros(n)
    exposure[rows], exposure_basis = native_exposure(obs.iloc[rows], total, exposure_key)
    if (exposure_basis == "gene_axis_row_sum" and measurement is not None
            and not measurement.measured_in_all(obs[measurement.sample_key].to_numpy()[admitted])[genes].all()):
        # Some admitted cells' samples lack genes of this axis, so their row sums are not shown to be
        # full native depth.
        exposure_basis = TRUNCATED
    reviewed_admitted = admitted.copy()
    admitted &= exposure >= min_umi
    lab = label.astype(str).to_numpy()
    names = sorted(set(lab[admitted]))
    if len(names) < 2:
        raise ValueError("A reference needs at least two admitted components")
    K = len(names)
    k_cell = np.where(admitted, pd.Categorical(lab, categories=names).codes, -1)
    member = k_cell >= 0

    # ---- component metadata
    meta = pd.DataFrame(index=pd.Index(names, name="name"))
    for field, key in (("parent", parent_key), ("compartment", compartment_key)):
        meta[field] = _per_component(obs[key].to_numpy()[member], k_cell[member], names, field) if key else None
    meta["review_status"], review_sources = component_review_status(obs, k_cell, K, level=annotation_level,
                                                                     review_status_key=review_status_key)
    if components is not None:
        declared = components.set_index("component")
        declared.index = declared.index.astype(str)
        unknown = sorted(set(declared.index) - set(names))
        if unknown:
            raise ValueError(f"Declared components have no admitted cells: {unknown[:5]}")
        for field in ("parent", "compartment"):
            if field in declared:
                values = declared[field].dropna()
                meta.loc[values.index, field] = values.astype(str).to_numpy()
        if "review_status" in declared:
            values = declared.review_status.dropna().astype(str)
            unknown = set(values) - set(REVIEW_STATUSES)
            if unknown:
                raise ValueError(f"review_status must be one of {REVIEW_STATUSES}; got {sorted(unknown)}")
            # A declaration adds caution; it cannot review what the annotation did not.
            meta.loc[values.index, "review_status"] = [
                REVIEW_STATUSES[max(REVIEW_STATUSES.index(a), REVIEW_STATUSES.index(b))]
                for a, b in zip(meta.loc[values.index, "review_status"], values)]

    # ---- per-component holdout by rule: min(max_cells, floor(fraction * n)) cells, split over donors
    # in proportion to their cells (largest remainder), so no donor is drained into the holdout;
    # fewer than holdout_min_cells is not held out (too few to test anything).
    donor = obs[donor_key].astype(str).to_numpy()
    holdout = np.zeros(n, bool)
    holdout_status = []
    for k in range(K):
        rows = np.flatnonzero(k_cell == k)
        n_hold = min(holdout_max_cells, int(holdout_fraction * len(rows)))
        if holdout_fraction <= 0 or holdout_max_cells <= 0:
            holdout_status.append("not_requested")
            continue
        if n_hold < holdout_min_cells:
            holdout_status.append("below_min_cells")
            continue
        donors, of = np.unique(donor[rows], return_inverse=True)
        share = n_hold * np.bincount(of) / len(rows)
        quota = np.floor(share).astype(int)
        order = rng.permutation(len(donors))  # random ties
        order = order[np.argsort((quota - share)[order], kind="stable")]
        quota[order[:n_hold - quota.sum()]] += 1
        picked = np.concatenate([rng.permutation(rows[of == d])[:quota[d]] for d in range(len(donors))])
        holdout[picked] = True
        holdout_status.append("sufficient" if len(set(donor[picked])) >= min_donors else "below_min_donors")
    fit = member & ~holdout

    # ---- one streamed pass
    groups = sorted(set(obs[group_key].astype(str).to_numpy()[fit])) if group_key else ["all"]
    g_cell = pd.Categorical(obs[group_key].astype(str), categories=groups).codes if group_key else np.zeros(n, int)
    combo = np.where(fit, k_cell * len(groups) + g_cell, -1).astype(np.int64)
    n_combo = K * len(groups)
    rows = np.flatnonzero(fit | holdout)
    jobs = [(path, layer, genes, block, combo, exposure, n_combo, holdout, chunk_rows)
            for block in np.array_split(rows, max(1, 2 * workers)) if len(block)]
    S = np.zeros((n_combo, len(genes)))
    held = []
    for sums, parts in _map_jobs(_stream, jobs, workers):
        S += sums
        held += parts
    n_combo_cells = np.bincount(combo[fit], minlength=n_combo).astype(np.float64)
    if measurement is not None:
        codes = measurement.sample_codes(obs[measurement.sample_key])
        N = np.zeros((len(measurement.samples), n_combo))
        np.add.at(N, (codes[fit], combo[fit]), 1)
        NM = N.T @ measurement.measured[:, genes].astype(np.float64)
    else:
        NM = np.repeat(n_combo_cells[:, None], len(genes), axis=1)

    # ---- measurement-aware profiles and fit axis
    Sk, NMk = (m.reshape(K, len(groups), -1).sum(1) for m in (S, NM))
    nk = n_combo_cells.reshape(K, len(groups)).sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        full = np.where(NMk > 0, Sk / NMk, np.nan).T  # genes x K
        frac = (NMk / nk[:, None]).T
    fit_gene = ((frac >= min_measured_fraction) & (NMk.T >= np.minimum(min_measured_cells, nk))).all(1)
    fit_gene &= np.isfinite(full).all(1) & (np.nansum(full, 1) > 0)
    if not fit_gene.any():
        raise ValueError("No gene is measured in enough training cells of every component")
    gene_names = var_names[genes].to_numpy()

    # ---- held-out validation genes: chosen on the fit axis by the RCTD DE rule, never fitted
    P = full[fit_gene]
    others = (P.sum(1, keepdims=True) - P) / max(K - 1, 1)
    lfc = np.log(P + 1e-9) - np.log(others + 1e-9)
    fully_measured = (frac[fit_gene] >= 0.95).all(1)
    heldout_rows, chosen = [], set()
    for k, name in enumerate(names):
        candidates = [i for i in np.flatnonzero((lfc[:, k] > RCTD_FC) & (P[:, k] > RCTD_EXPR) & fully_measured)
                      if i not in chosen]
        n_take = min(heldout_genes_per_component, len(candidates) // 5)  # keep >= 80% of a type's DE genes
        for i in (rng.choice(candidates, n_take, replace=False) if n_take else []):
            chosen.add(int(i))
            heldout_rows.append({"gene": gene_names[fit_gene][i], "selection": "type_marker", "source_component": name,
                                 "logFC_vs_other_mean": float(lfc[i, k]), "profile": float(P[i, k])})
    background = np.flatnonzero(~((lfc > RCTD_FC) & (P > RCTD_EXPR)).any(1) & fully_measured & (P.mean(1) > 1e-4))
    for i in rng.choice(background, min(heldout_random_genes, len(background)), replace=False):
        heldout_rows.append({"gene": gene_names[fit_gene][i], "selection": "random_expressed_nonDE", "source_component": "",
                             "logFC_vs_other_mean": float(lfc[i].max()), "profile": float(P[i].mean())})
    heldout_genes = pd.DataFrame(heldout_rows, columns=["gene", "selection", "source_component", "logFC_vs_other_mean", "profile"])

    # ---- outputs
    profile_path = out / "reference_profiles.npz"
    np.savez(profile_path, profiles=P, gene_names=gene_names[fit_gene].astype(str), cell_type_names=np.asarray(names, dtype=str))
    np.savez_compressed(out / "component_sums.npz", sum_x_over_numi=S, measured_cells=NM, n_cells=n_combo_cells,
                        components=np.asarray(names, dtype=str), groups=np.asarray(groups, dtype=str),
                        gene_names=gene_names.astype(str), in_fit_axis=fit_gene)
    heldout_genes.to_csv(out / "heldout_validation_genes.tsv", sep="\t", index=False)
    import anndata as ad

    hold_rows = np.concatenate([r for r, _ in held]) if held else np.empty(0, int)
    order = np.argsort(hold_rows)
    hold_x = sparse.vstack([x for _, x in held]).tocsr()[order] if held else sparse.csr_matrix((0, len(genes)))
    hold_rows = hold_rows[order]
    hold_var = pd.DataFrame(index=pd.Index(gene_names.astype(str)))
    if measurement is not None:
        # Pseudo-spots pool holdout cells across samples; only these genes are measured in all of them.
        hold_var["measured_in_all_holdout_samples"] = measurement.measured_in_all(
            obs[measurement.sample_key].to_numpy()[hold_rows])[genes]
    ad.AnnData(hold_x, obs=pd.DataFrame({"component": np.asarray(names)[k_cell[hold_rows]], "donor": donor[hold_rows],
                                         "native_full_exposure": exposure[hold_rows]}, index=obs.index[hold_rows].astype(str)),
               var=hold_var).write_h5ad(out / "holdout_cells.h5ad")
    coverage = pd.DataFrame({"label": lab, "n_cells": 1, "n_excluded_technical": technical,
                             "n_not_admitted": ~technical & ~reviewed_admitted,
                             "n_below_min_umi": reviewed_admitted & ~admitted, "n_training": fit, "n_holdout": holdout}
                            ).groupby("label").sum().astype(int)
    coverage["decision"] = np.select([coverage.n_training + coverage.n_holdout > 0, coverage.n_excluded_technical == coverage.n_cells],
                                     ["component", "excluded_technical"], "not_admitted")
    coverage.reset_index().to_csv(out / "reference_coverage.tsv", sep="\t", index=False)

    chemistry = obs[chemistry_key].astype(str).to_numpy() if chemistry_key else np.full(n, "undeclared", dtype=object)
    mean_exposure = pd.Series(exposure[fit]).groupby([k_cell[fit], chemistry[fit]]).mean()
    native_preserved = exposure_basis != TRUNCATED
    train_donors = pd.Series(donor[fit]).groupby(k_cell[fit]).agg(set)
    hold_donors = pd.Series(donor[holdout]).groupby(k_cell[holdout]).agg(set)
    hold_cells = np.bincount(k_cell[holdout], minlength=K)
    entries = []
    for k, name in enumerate(names):
        trained, held_by = train_donors.get(k, set()), hold_donors.get(k, set())
        entries.append({"name": name, "parent": meta.parent.iloc[k], "compartment": meta.compartment.iloc[k],
                        "review_status": meta.review_status.iloc[k], "training_cells": int(nk[k]),
                        "training_donors": len(trained), "holdout_cells": int(hold_cells[k]),
                        "holdout_donors": len(held_by), "holdout_donors_in_training": len(held_by & trained),
                        "holdout_status": holdout_status[k],
                        "mean_native_exposure_by_chemistry": ({c: float(v) for c, v in mean_exposure.loc[k].items()}
                                                              if native_preserved else None),
                        "measured_fraction_min": float(frac[fit_gene, k].min())})
    gaps = []
    if not native_preserved:
        gaps.append("native exposure not established: admitted cells' samples lack genes of the axis, so their "
                    "row sums may miss native depth; declare --exposure-key. RCTD fits refuse this package")
    for c in entries:
        if c["review_status"] == "evidence_pending":
            gaps.append(f"{c['name']}: evidence_pending component; production fits refuse it")
        elif c["review_status"] != "reviewed":
            gaps.append(f"{c['name']}: {c['review_status']} component; its fine abundance is not interpretable")
        if c["training_cells"] < min_cells:
            gaps.append(f"{c['name']}: {c['training_cells']} training cells (< {min_cells})")
        if c["training_donors"] < min_donors:
            gaps.append(f"{c['name']}: {c['training_donors']} training donor(s); donor-specific profile")
        if c["holdout_status"] in {"below_min_cells", "below_min_donors"}:
            gaps.append(f"{c['name']}: holdout {c['holdout_status']}; pseudo-spot identifiability "
                        + ("not assessable" if c["holdout_status"] == "below_min_cells" else "donor-confounded"))
    gaps.append("whole-tissue reference coverage is not assessed by the builder")
    manifest = {
        "schema_version": SCHEMA, "unit": UNIT, "native_full_exposure_preserved": native_preserved,
        "exposure_basis": exposure_basis,
        "profiles_path": profile_path.name, "profiles_sha256": sha256_file(profile_path),
        "estimator": "per gene: sum over training cells of type whose sample measured the gene of x/nUMI, "
                     "divided by their number; nUMI = full native exposure",
        "components": entries, "coverage_complete": False, "gaps": gaps,
        "holdout_rule": {"cells_per_component": "min(max_cells, floor(fraction * n_admitted))",
                         "fraction": holdout_fraction, "max_cells": holdout_max_cells, "min_cells": holdout_min_cells,
                         "min_donors": min_donors, "seed": seed,
                         "assignment": "per donor, in proportion to its cells (largest remainder)",
                         "scope": "cell holdout within donors: holdout donors also train (holdout_donors_in_training)",
                         "below_min_cells": "not held out; all cells train", "path": "holdout_cells.h5ad",
                         "n_cells": int(holdout.sum())},
        "heldout_validation_genes": {"path": "heldout_validation_genes.tsv", "n": len(heldout_genes),
                                     "use": "pass as spatial.rctd.exclude_genes; never fitted"},
        "measurement": {"declared": measurement is not None,
                        "min_measured_fraction_all_components": min_measured_fraction,
                        "min_measured_cells_all_components": min_measured_cells,
                        "fit_axis_genes": int(fit_gene.sum()), "gene_axis_genes": len(genes)},
        "admission": {"label_key": label_key, "status_key": status_key, "admit_status": list(admit_status),
                      "never_admitted": "excluded_technical", "reviewed_key": reviewed_key if reviewed_key in obs else None,
                      "technical_key": technical_key, "require_true": list(require_true), "min_umi": min_umi,
                      "exposure": exposure_key or "row sum of gene-expression features",
                      "review_status_key": review_status_key, "annotation_level": annotation_level,
                      "review_status_sources": review_sources, "donor_key": donor_key, "group_key": group_key,
                      "chemistry_key": chemistry_key, "n_admitted_cells": int(admitted.sum()),
                      "n_training_cells": int(fit.sum()), "n_excluded_technical": int(technical.sum())},
        "cap_policy": "no per-type cap: every admitted non-holdout cell enters the pooled mean",
        "component_sums": {"path": "component_sums.npz", "groups": group_key or "all",
                           "use": "exact re-aggregation of merged components (same cells, same estimator)"},
        "inputs": {"counts_h5ad": {"path": str(path.resolve()), "size_bytes": path.stat().st_size,
                                   "var_sha256": var_sha256(var_names), "layer": layer or "X"}},
        "producer": {"scoop_version": __version__, "module_sha256": sha256_file(__file__)},
        "elapsed_seconds": perf_counter() - started,
    }
    if measurement is not None:
        manifest["measurement"].update({"sample_key": measurement.sample_key, "var_sha256": measurement.var_sha256})
    atomic_json(out / "reference_manifest.json", manifest)
    return manifest
