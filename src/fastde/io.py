from __future__ import annotations

from pathlib import Path
import json

import pandas as pd


def read_pseudobulk_dir(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = Path(path)
    counts = pd.read_csv(path / "counts.tsv", sep="\t")
    meta = pd.read_csv(path / "metadata.tsv", sep="\t")
    if "pseudobulk_id" not in counts.columns or "pseudobulk_id" not in meta.columns:
        raise ValueError("counts.tsv and metadata.tsv must contain pseudobulk_id")
    counts = counts.set_index("pseudobulk_id")
    meta = meta.set_index("pseudobulk_id")
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    return counts, meta


def write_json(obj, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path
