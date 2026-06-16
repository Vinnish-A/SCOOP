from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from .abundance_data import AbundanceTable
from .abundance_design import build_covariate_matrix, parse_covariates
from .abundance_loss import cox_partial_nll
from .abundance_metrics import benjamini_hochberg, binary_metrics, concordance_index, correlation_z_p, multiclass_metrics, welch_z_p
from .abundance_result import AbundanceResult


@dataclass
class BagDataset:
    sample_ids: list[str]
    bags: list[np.ndarray]
    celltype_names: list[str]
    covariates: np.ndarray | None = None


class GatedAttention(nn.Module):
    """Gated attention pooling matching scSurvival's bag-level design."""

    def __init__(self, hidden_dim: int):
        super().__init__()
        self.v = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.u = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.w = nn.Linear(hidden_dim, 1, bias=False)

    def forward(self, h: torch.Tensor) -> torch.Tensor:
        return self.w(torch.tanh(self.v(h)) * torch.sigmoid(self.u(h)))


class ScSurvivalMILNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dim: int = 64,
        dropout: float = 0.1,
        covariate_dim: int = 0,
    ):
        super().__init__()
        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )
        self.attention = GatedAttention(hidden_dim)
        head_in = hidden_dim + int(covariate_dim)
        self.head = nn.Sequential(
            nn.Linear(head_in, max(1, head_in // 2)),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(max(1, head_in // 2), output_dim),
        )
        self.instance_head = nn.Linear(hidden_dim, output_dim)

    def encode_instances(self, x: torch.Tensor) -> torch.Tensor:
        return self.feature_extractor(x)

    def forward_bag(self, x: torch.Tensor, covariates: torch.Tensor | None = None) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self.encode_instances(x)
        attention = torch.softmax(self.attention(h), dim=0)
        pooled = torch.sum(h * attention, dim=0)
        if covariates is not None:
            pooled = torch.cat([pooled, covariates], dim=0)
        out = self.head(pooled)
        instance_out = self.instance_head(h)
        return out, attention.view(-1), instance_out

    def score_prototypes(self, prototypes: torch.Tensor) -> torch.Tensor:
        h = self.encode_instances(prototypes)
        return self.instance_head(h)


def _cox_loss_torch(risk: torch.Tensor, time: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    order = torch.argsort(time.view(-1), descending=True)
    risk = risk.view(-1)[order]
    event = event.view(-1)[order]
    log_risk = torch.logcumsumexp(risk, dim=0)
    observed = event > 0
    if not torch.any(observed):
        raise ValueError("Cox partial likelihood requires at least one event")
    return -(risk[observed] - log_risk[observed]).mean()


def _cox_rank_loss_torch(risk: torch.Tensor, time: torch.Tensor, event: torch.Tensor, margin: float = 1.0) -> torch.Tensor:
    risk = risk.view(-1)
    time = time.view(-1)
    event = event.view(-1)
    losses = []
    for i in range(len(time)):
        if event[i] <= 0:
            continue
        comparable = time[i] < time
        if torch.any(comparable):
            losses.append(F.relu(margin - (risk[i] - risk[comparable])).mean())
    if not losses:
        return risk.sum() * 0.0
    return torch.stack(losses).mean()


def _prepare_covariates(table: AbundanceTable, covariates: str | list[str] | None) -> pd.DataFrame:
    return build_covariate_matrix(table.metadata, parse_covariates(covariates))


def build_bags_from_counts(table: AbundanceTable, *, max_instances_per_sample: int = 2000, random_seed: int = 0) -> BagDataset:
    rng = np.random.default_rng(random_seed)
    celltypes = list(table.counts.columns)
    eye = np.eye(len(celltypes), dtype=np.float32)
    bags: list[np.ndarray] = []
    for _, row in table.counts.iterrows():
        counts = row.to_numpy(dtype=int)
        total = int(counts.sum())
        if total <= 0:
            raise ValueError("sample has zero cells")
        if total > max_instances_per_sample:
            probs = counts / total
            counts = rng.multinomial(max_instances_per_sample, probs)
        bags.append(np.repeat(eye, counts, axis=0))
    return BagDataset(sample_ids=list(table.counts.index.astype(str)), bags=bags, celltype_names=celltypes)


def build_bags_from_h5ad(
    input_h5ad: str | Path,
    table: AbundanceTable,
    sample_key: str,
    celltype_key: str,
    *,
    max_instances_per_sample: int = 2000,
    random_seed: int = 0,
) -> BagDataset:
    rng = np.random.default_rng(random_seed)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    obs = adata.obs[[sample_key, celltype_key]].dropna().copy()
    obs[sample_key] = obs[sample_key].astype(str)
    obs[celltype_key] = obs[celltype_key].astype(str)
    celltypes = list(table.counts.columns.astype(str))
    mapping = {name: i for i, name in enumerate(celltypes)}
    eye = np.eye(len(celltypes), dtype=np.float32)
    bags: list[np.ndarray] = []
    for sample in table.counts.index.astype(str):
        sample_obs = obs[(obs[sample_key] == sample) & obs[celltype_key].isin(mapping)]
        idx = sample_obs[celltype_key].map(mapping).to_numpy(dtype=int)
        if len(idx) == 0:
            raise ValueError(f"sample {sample!r} has no usable cell-type labels")
        if len(idx) > max_instances_per_sample:
            idx = rng.choice(idx, size=max_instances_per_sample, replace=False)
        bags.append(eye[idx])
    if getattr(adata, "isbacked", False):
        adata.file.close()
    return BagDataset(sample_ids=list(table.counts.index.astype(str)), bags=bags, celltype_names=celltypes)


def _fit_mil(
    dataset: BagDataset,
    *,
    mode: str,
    y: np.ndarray | None = None,
    time: np.ndarray | None = None,
    event: np.ndarray | None = None,
    n_classes: int | None = None,
    covariates: np.ndarray | None = None,
    hidden_dim: int = 64,
    dropout: float = 0.1,
    learning_rate: float = 1e-3,
    weight_decay: float = 1e-4,
    max_epochs: int = 500,
    random_seed: int = 0,
    survival_loss: str = "cox",
) -> tuple[ScSurvivalMILNet, np.ndarray, list[dict[str, float]]]:
    torch.manual_seed(random_seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    output_dim = int(n_classes or 1)
    covariate_dim = 0 if covariates is None else int(covariates.shape[1])
    model = ScSurvivalMILNet(len(dataset.celltype_names), output_dim, hidden_dim=hidden_dim, dropout=dropout, covariate_dim=covariate_dim).to(device)
    bags = [torch.tensor(bag, dtype=torch.float32, device=device) for bag in dataset.bags]
    cov_tensors = None if covariates is None else [torch.tensor(row, dtype=torch.float32, device=device) for row in covariates]
    opt = torch.optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    history: list[dict[str, float]] = []
    y_t = None if y is None else torch.tensor(y, dtype=torch.long if mode == "multiclass" else torch.float32, device=device)
    time_t = None if time is None else torch.tensor(time, dtype=torch.float32, device=device)
    event_t = None if event is None else torch.tensor(event, dtype=torch.float32, device=device)
    for epoch in range(int(max_epochs)):
        model.train()
        opt.zero_grad()
        outputs = []
        for i, bag in enumerate(bags):
            cov = None if cov_tensors is None else cov_tensors[i]
            out, _, _ = model.forward_bag(bag, cov)
            outputs.append(out)
        logits = torch.stack(outputs)
        if mode == "survival":
            risk = logits.view(-1)
            assert time_t is not None and event_t is not None
            if survival_loss == "cox_rank":
                loss = _cox_rank_loss_torch(risk, time_t, event_t)
            elif survival_loss == "cox_plus_rank":
                loss = _cox_loss_torch(risk, time_t, event_t) + 0.1 * _cox_rank_loss_torch(risk, time_t, event_t)
            else:
                loss = _cox_loss_torch(risk, time_t, event_t)
        elif mode in {"binary", "condition"}:
            assert y_t is not None
            loss = F.binary_cross_entropy_with_logits(logits.view(-1), y_t.float())
        elif mode == "multiclass":
            assert y_t is not None
            loss = F.cross_entropy(logits, y_t.long())
        elif mode == "continuous":
            assert y_t is not None
            loss = F.mse_loss(logits.view(-1), y_t.float())
        else:
            raise ValueError(f"unknown MIL mode: {mode}")
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        if epoch % 25 == 0 or epoch == max_epochs - 1:
            history.append({"epoch": float(epoch), "loss": float(loss.detach().cpu())})
    model.eval()
    with torch.no_grad():
        pred = []
        for i, bag in enumerate(bags):
            cov = None if cov_tensors is None else cov_tensors[i]
            out, _, _ = model.forward_bag(bag, cov)
            pred.append(out.detach().cpu().numpy())
    return model, np.vstack(pred), history


def _prototype_scores(model: ScSurvivalMILNet, n_celltypes: int) -> np.ndarray:
    device = next(model.parameters()).device
    with torch.no_grad():
        eye = torch.eye(n_celltypes, dtype=torch.float32, device=device)
        scores = model.score_prototypes(eye).detach().cpu().numpy()
    return scores


def _mil_binary_results(table: AbundanceTable, y: np.ndarray, scores: np.ndarray, positive_label: str, negative_label: str, mode: str) -> pd.DataFrame:
    rows = []
    props = table.proportions
    for i, cell_type in enumerate(props.columns):
        pos = props.iloc[y == 1, i].to_numpy(dtype=float)
        neg = props.iloc[y == 0, i].to_numpy(dtype=float)
        z, pvalue = welch_z_p(pos, neg)
        c = float(scores[i])
        rows.append(
            {
                "cell_type": cell_type,
                "coef_or_importance": c,
                "direction": "positive" if c >= 0 else "negative",
                "odds_ratio_approx": float(np.exp(np.clip(c, -20, 20))),
                "mean_prop_positive": float(pos.mean()) if len(pos) else np.nan,
                "mean_prop_negative": float(neg.mean()) if len(neg) else np.nan,
                "delta_prop": float(pos.mean() - neg.mean()) if len(pos) and len(neg) else np.nan,
                "z": z,
                "pvalue": pvalue,
                "model": "scsurvival_mil_attention",
                "mode": mode,
                "positive_label": positive_label,
                "negative_label": negative_label,
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].to_numpy(dtype=float))
    return out.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)


def _mil_survival_results(table: AbundanceTable, time: np.ndarray, event: np.ndarray, risk: np.ndarray, scores: np.ndarray) -> pd.DataFrame:
    rows = []
    props = table.proportions
    high = risk >= np.median(risk)
    cindex = concordance_index(time, risk, event)
    for i, cell_type in enumerate(props.columns):
        feature = props.iloc[:, i].to_numpy(dtype=float)
        z, pvalue = correlation_z_p(feature, risk)
        c = float(scores[i])
        rows.append(
            {
                "cell_type": cell_type,
                "coef_or_importance": c,
                "hazard_ratio_approx": float(np.exp(np.clip(c, -20, 20))),
                "direction": "high_risk" if c >= 0 else "low_risk",
                "mean_prop_high_risk": float(feature[high].mean()) if high.any() else np.nan,
                "mean_prop_low_risk": float(feature[~high].mean()) if (~high).any() else np.nan,
                "z": z,
                "pvalue": pvalue,
                "concordance_index": cindex,
                "n_events": int(event.sum()),
                "n_samples": int(len(time)),
                "model": "scsurvival_mil_attention",
                "mode": "survival",
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].to_numpy(dtype=float))
    return out.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)


def run_mil_abundance(
    *,
    mode: str,
    table: AbundanceTable,
    output_dir: str | Path,
    input_h5ad: str | Path | None = None,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type_lvl3",
    label_col: str | None = None,
    positive_label: str | None = None,
    negative_label: str | None = None,
    reference_level: str | None = None,
    time_col: str | None = None,
    event_col: str | None = None,
    value_col: str | None = None,
    covariates: str | list[str] | None = None,
    max_instances_per_sample: int = 2000,
    hidden_dim: int = 64,
    dropout: float = 0.1,
    learning_rate: float = 1e-3,
    weight_decay: float = 1e-4,
    max_epochs: int = 500,
    random_seed: int = 0,
    survival_loss: str = "cox",
    manifest_factory: Any | None = None,
) -> AbundanceResult:
    from .abundance_design import encode_binary_labels, encode_multiclass_labels, encode_survival

    if manifest_factory is None:
        raise ValueError("manifest_factory is required")
    dataset = (
        build_bags_from_h5ad(input_h5ad, table, sample_key, celltype_key, max_instances_per_sample=max_instances_per_sample, random_seed=random_seed)
        if input_h5ad is not None
        else build_bags_from_counts(table, max_instances_per_sample=max_instances_per_sample, random_seed=random_seed)
    )
    cov_df = _prepare_covariates(table, covariates)
    cov = cov_df.to_numpy(dtype=np.float32) if cov_df.shape[1] else None
    if mode == "binary":
        if not label_col or positive_label is None or negative_label is None:
            raise ValueError("binary mode requires label_col, positive_label, and negative_label")
        y, keep = encode_binary_labels(table.metadata, label_col, positive_label, negative_label)
        keep_pos = [table.counts.index.get_loc(s) for s in keep]
        table = AbundanceTable(table.counts.loc[keep], table.proportions.loc[keep], table.metadata.loc[keep], table.sample_key, table.celltype_key)
        dataset = BagDataset([dataset.sample_ids[i] for i in keep_pos], [dataset.bags[i] for i in keep_pos], dataset.celltype_names)
        cov = None if cov is None else cov[keep_pos]
        model, logits, history = _fit_mil(dataset, mode="binary", y=y, covariates=cov, hidden_dim=hidden_dim, dropout=dropout, learning_rate=learning_rate, weight_decay=weight_decay, max_epochs=max_epochs, random_seed=random_seed)
        logit = logits.reshape(-1)
        metrics = binary_metrics(y, logit)
        proto = _prototype_scores(model, table.n_celltypes).reshape(-1)
        results = _mil_binary_results(table, y.astype(int), proto, str(positive_label), str(negative_label), "binary")
        predictions = pd.DataFrame({"sample_id": table.counts.index, "label": y.astype(int), "logit": logit, "probability": 1 / (1 + np.exp(-logit))})
        manifest = manifest_factory(mode="binary", table=table, sample_key=sample_key, celltype_key=celltype_key, transform="bag_onehot", n_classes=2, metrics=metrics, model={"backend": "scsurvival_mil", "history": history, "survival_loss": None})
    elif mode == "survival":
        if not time_col or not event_col:
            raise ValueError("survival mode requires time_col and event_col")
        time, event, keep = encode_survival(table.metadata, time_col, event_col)
        keep_pos = [table.counts.index.get_loc(s) for s in keep]
        table = AbundanceTable(table.counts.loc[keep], table.proportions.loc[keep], table.metadata.loc[keep], table.sample_key, table.celltype_key)
        dataset = BagDataset([dataset.sample_ids[i] for i in keep_pos], [dataset.bags[i] for i in keep_pos], dataset.celltype_names)
        cov = None if cov is None else cov[keep_pos]
        model, logits, history = _fit_mil(dataset, mode="survival", time=time, event=event, covariates=cov, hidden_dim=hidden_dim, dropout=dropout, learning_rate=learning_rate, weight_decay=weight_decay, max_epochs=max_epochs, random_seed=random_seed, survival_loss=survival_loss)
        risk = logits.reshape(-1)
        metrics = {"cox_partial_nll": cox_partial_nll(risk, time, event), "concordance_index": concordance_index(time, risk, event)}
        proto = _prototype_scores(model, table.n_celltypes).reshape(-1)
        results = _mil_survival_results(table, time, event, risk, proto)
        predictions = pd.DataFrame({"sample_id": table.counts.index, "time": time, "event": event.astype(int), "risk": risk})
        manifest = manifest_factory(mode="survival", table=table, sample_key=sample_key, celltype_key=celltype_key, transform="bag_onehot", n_events=int(event.sum()), metrics=metrics, model={"backend": "scsurvival_mil", "history": history, "survival_loss": survival_loss})
    elif mode == "multiclass":
        if not label_col:
            raise ValueError("multiclass mode requires label_col")
        y, classes, keep, reference_idx = encode_multiclass_labels(table.metadata, label_col, reference_level)
        keep_pos = [table.counts.index.get_loc(s) for s in keep]
        table = AbundanceTable(table.counts.loc[keep], table.proportions.loc[keep], table.metadata.loc[keep], table.sample_key, table.celltype_key)
        dataset = BagDataset([dataset.sample_ids[i] for i in keep_pos], [dataset.bags[i] for i in keep_pos], dataset.celltype_names)
        cov = None if cov is None else cov[keep_pos]
        model, logits, history = _fit_mil(dataset, mode="multiclass", y=y, n_classes=len(classes), covariates=cov, hidden_dim=hidden_dim, dropout=dropout, learning_rate=learning_rate, weight_decay=weight_decay, max_epochs=max_epochs, random_seed=random_seed)
        metrics = multiclass_metrics(y, logits, classes)
        proto = _prototype_scores(model, table.n_celltypes)
        if reference_idx is None:
            reference_idx = 0
        rows = []
        for class_idx, cls in enumerate(classes):
            if class_idx == reference_idx:
                continue
            contrast = proto[:, class_idx] - proto[:, reference_idx]
            for i, cell_type in enumerate(table.proportions.columns):
                cur = table.proportions.iloc[y == class_idx, i].to_numpy(dtype=float)
                ref = table.proportions.iloc[y == reference_idx, i].to_numpy(dtype=float)
                z, pvalue = welch_z_p(cur, ref)
                rows.append({"cell_type": cell_type, "class_or_contrast": f"{cls}_vs_{classes[reference_idx]}", "coef_or_importance": float(contrast[i]), "direction": cls if contrast[i] >= 0 else classes[reference_idx], "mean_prop_class": float(cur.mean()) if len(cur) else np.nan, "mean_prop_reference": float(ref.mean()) if len(ref) else np.nan, "delta_prop": float(cur.mean() - ref.mean()) if len(cur) and len(ref) else np.nan, "z": z, "pvalue": pvalue, "model": "scsurvival_mil_attention", "mode": "multiclass"})
        results = pd.DataFrame(rows)
        results["padj"] = benjamini_hochberg(results["pvalue"].to_numpy(dtype=float))
        results = results.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)
        predictions = pd.DataFrame({"sample_id": table.counts.index, "label": [classes[int(i)] for i in y], "predicted_label": [classes[int(i)] for i in logits.argmax(axis=1)]})
        probs = np.exp(logits - logits.max(axis=1, keepdims=True))
        probs = probs / probs.sum(axis=1, keepdims=True)
        for i, cls in enumerate(classes):
            predictions[f"prob_{cls}"] = probs[:, i]
        manifest = manifest_factory(mode="multiclass", table=table, sample_key=sample_key, celltype_key=celltype_key, transform="bag_onehot", n_classes=len(classes), metrics=metrics, model={"backend": "scsurvival_mil", "history": history, "survival_loss": None})
    elif mode == "continuous":
        if not value_col or value_col not in table.metadata.columns:
            raise ValueError("continuous mode requires value_col in metadata")
        values = pd.to_numeric(table.metadata[value_col], errors="coerce")
        keep = table.metadata.index[values.notna()]
        keep_pos = [table.counts.index.get_loc(s) for s in keep]
        y = values.loc[keep].to_numpy(dtype=float)
        table = AbundanceTable(table.counts.loc[keep], table.proportions.loc[keep], table.metadata.loc[keep], table.sample_key, table.celltype_key)
        dataset = BagDataset([dataset.sample_ids[i] for i in keep_pos], [dataset.bags[i] for i in keep_pos], dataset.celltype_names)
        cov = None if cov is None else cov[keep_pos]
        model, logits, history = _fit_mil(dataset, mode="continuous", y=y, covariates=cov, hidden_dim=hidden_dim, dropout=dropout, learning_rate=learning_rate, weight_decay=weight_decay, max_epochs=max_epochs, random_seed=random_seed)
        pred = logits.reshape(-1)
        z, pvalue = correlation_z_p(pred, y)
        metrics = {"prediction_correlation_z": z, "prediction_correlation_pvalue": pvalue}
        proto = _prototype_scores(model, table.n_celltypes).reshape(-1)
        rows = []
        for i, cell_type in enumerate(table.proportions.columns):
            fz, fp = correlation_z_p(table.proportions.iloc[:, i].to_numpy(dtype=float), y)
            rows.append({"cell_type": cell_type, "coef_or_importance": float(proto[i]), "direction": "positive" if proto[i] >= 0 else "negative", "z": fz, "pvalue": fp, "model": "scsurvival_mil_attention", "mode": "continuous"})
        results = pd.DataFrame(rows)
        results["padj"] = benjamini_hochberg(results["pvalue"].to_numpy(dtype=float))
        predictions = pd.DataFrame({"sample_id": table.counts.index, "value": y, "prediction": pred})
        manifest = manifest_factory(mode="continuous", table=table, sample_key=sample_key, celltype_key=celltype_key, transform="bag_onehot", metrics=metrics, model={"backend": "scsurvival_mil", "history": history, "survival_loss": None})
    else:
        raise ValueError(f"unknown MIL abundance mode: {mode}")
    result = AbundanceResult(mode=mode, results=results, predictions=predictions, metrics=metrics, manifest=manifest)
    result.write(output_dir)
    return result
