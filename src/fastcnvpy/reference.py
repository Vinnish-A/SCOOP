from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_gene_metadata(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        frame = pd.read_csv(path, sep="\t", low_memory=False)
    elif suffix == ".csv":
        frame = pd.read_csv(path)
    elif suffix in {".parquet", ".pq"}:
        frame = pd.read_parquet(path)
    else:
        raise ValueError(f"unsupported gene metadata format: {path}")
    return normalize_gene_metadata(frame)


def normalize_gene_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    normalized_required = {"hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm"}
    if normalized_required.issubset(frame.columns):
        out = frame.copy()
        out["chromosome_num"] = pd.to_numeric(out["chromosome_num"], errors="coerce")
        out["start_position"] = pd.to_numeric(out["start_position"], errors="coerce")
        out["end_position"] = pd.to_numeric(out["end_position"], errors="coerce")
        out = out.dropna(subset=["hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm"])
        out["chromosome_num"] = out["chromosome_num"].astype(int)
        out["hgnc_symbol"] = out["hgnc_symbol"].astype(str)
        out["chr_arm"] = _arm_suffix(out["chr_arm"], out["chromosome_num"])
        return out[["hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm"]].reset_index(drop=True)

    required = {
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position",
        "gene_biotype",
        "chr_arm",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"gene metadata missing columns: {sorted(missing)}")

    out = frame.copy()
    out["chromosome_name"] = out["chromosome_name"].astype(str)
    keep = (
        out["gene_biotype"].isin(["protein_coding", "lncRNA"])
        & out["chromosome_name"].isin([str(i) for i in range(1, 23)] + ["X"])
        & out["hgnc_symbol"].fillna("").astype(str).ne("")
    )
    out = out.loc[keep].copy()
    out["chromosome_num"] = out["chromosome_name"].replace({"X": "23"}).astype(float).astype(int)
    out["start_position"] = pd.to_numeric(out["start_position"], errors="coerce")
    out["end_position"] = pd.to_numeric(out["end_position"], errors="coerce")
    out = out.dropna(subset=["start_position", "end_position", "chr_arm"])
    out["hgnc_symbol"] = out["hgnc_symbol"].astype(str)
    out["chr_arm"] = _arm_suffix(out["chr_arm"], out["chromosome_num"])
    out = out.drop_duplicates(["hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm"])
    return out[["hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm"]].reset_index(drop=True)


def _arm_suffix(arms: pd.Series, chromosomes: pd.Series) -> pd.Series:
    values = arms.astype(str).str.strip()
    suffix = values.str.extract(r"^(?:chr)?(?:(?:[0-9]+|X)\.?)?([pq])$", expand=False)
    prefix = values.str.extract(r"^(?:chr)?([0-9]+|X)\.?[pq]$", expand=False).where(lambda x: x.ne("X"), "23")
    mismatch = prefix.notna() & pd.to_numeric(prefix, errors="coerce").ne(chromosomes)
    if suffix.isna().any() or mismatch.any():
        raise ValueError("chr_arm must be p/q or a matching chromosome-qualified arm (e.g. 19.q)")
    return suffix


def chromosome_arm_order() -> list[str]:
    return [f"{chrom}.{arm}" for chrom in range(1, 23) for arm in ("p", "q")] + ["X.p", "X.q"]


def window_chrom_arm(window_name: str) -> str:
    chrom, rest = window_name.split(".", 1)
    arm = rest[0]
    if chrom == "23":
        chrom = "X"
    return f"{chrom}.{arm}"
