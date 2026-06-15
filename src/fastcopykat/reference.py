from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


_GENE_ALIASES = {
    "gene": ("gene", "hgnc_symbol", "symbol", "Gene", "GeneSymbol", "ensembl_gene_id"),
    "chrom": ("chrom", "chromosome", "chromosome_name", "chr"),
    "start": ("start", "start_position", "chrompos", "position"),
    "end": ("end", "end_position"),
    "abspos": ("abspos", "abs_pos"),
}


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    raise ValueError(f"unsupported table format: {path}")


def _pick_column(columns: pd.Index, aliases: tuple[str, ...], *, required: bool = True) -> str | None:
    for alias in aliases:
        if alias in columns:
            return alias
    if required:
        raise ValueError(f"missing one of columns: {aliases}")
    return None


def normalize_gene_annotation(annotation: pd.DataFrame) -> pd.DataFrame:
    """Return CopyKAT-style gene coordinates with stable column names."""

    gene_col = _pick_column(annotation.columns, _GENE_ALIASES["gene"])
    chrom_col = _pick_column(annotation.columns, _GENE_ALIASES["chrom"])
    start_col = _pick_column(annotation.columns, _GENE_ALIASES["start"])
    end_col = _pick_column(annotation.columns, _GENE_ALIASES["end"], required=False)
    abspos_col = _pick_column(annotation.columns, _GENE_ALIASES["abspos"], required=False)

    out = pd.DataFrame(
        {
            "gene": annotation[gene_col].astype(str),
            "chrom": annotation[chrom_col].astype(str).str.replace("^chr", "", regex=True),
            "start": pd.to_numeric(annotation[start_col], errors="coerce"),
        }
    )
    if end_col is None:
        out["end"] = out["start"]
    else:
        out["end"] = pd.to_numeric(annotation[end_col], errors="coerce")
    if abspos_col is None:
        out["abspos"] = np.nan
    else:
        out["abspos"] = pd.to_numeric(annotation[abspos_col], errors="coerce")
    out["chrompos"] = ((out["start"] + out["end"]) / 2).round().astype("float64")
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["gene", "chrom", "start", "end", "chrompos"])
    out = out[out["gene"].ne("") & out["chrom"].ne("")]
    out = out.drop_duplicates("gene", keep="first")

    chrom_rank = _chromosome_rank(out["chrom"])
    chrom_offsets = _chromosome_offsets(out, chrom_rank)
    missing_abspos = out["abspos"].isna()
    out.loc[missing_abspos, "abspos"] = (
        out.loc[missing_abspos, "chrom"].map(chrom_offsets).astype(float) + out.loc[missing_abspos, "chrompos"]
    )
    out["_chrom_rank"] = chrom_rank
    out = out.sort_values(["_chrom_rank", "chrompos", "gene"]).drop(columns=["_chrom_rank"]).reset_index(drop=True)
    return out[["gene", "chrom", "start", "end", "chrompos", "abspos"]]


def normalize_bins(bins: pd.DataFrame) -> pd.DataFrame:
    chrom_col = _pick_column(bins.columns, ("chrom", "chromosome", "chromosome_name", "chr"))
    chrompos_col = _pick_column(bins.columns, ("chrompos", "start", "position"))
    abspos_col = _pick_column(bins.columns, ("abspos", "abs_pos"), required=False)
    out = pd.DataFrame(
        {
            "chrom": bins[chrom_col].astype(str).str.replace("^chr", "", regex=True),
            "chrompos": pd.to_numeric(bins[chrompos_col], errors="coerce"),
        }
    )
    if abspos_col is None:
        out["abspos"] = np.nan
    else:
        out["abspos"] = pd.to_numeric(bins[abspos_col], errors="coerce")
    out = out.dropna(subset=["chrom", "chrompos"])
    chrom_rank = _chromosome_rank(out["chrom"])
    if out["abspos"].isna().any():
        offsets = _chromosome_offsets_from_positions(out, chrom_rank)
        missing = out["abspos"].isna()
        out.loc[missing, "abspos"] = out.loc[missing, "chrom"].map(offsets).astype(float) + out.loc[missing, "chrompos"]
    out["_chrom_rank"] = chrom_rank
    return out.sort_values(["_chrom_rank", "chrompos"]).drop(columns=["_chrom_rank"]).reset_index(drop=True)


def build_bins_from_annotation(annotation: pd.DataFrame, bin_size: int = 220_000) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    chrom_rank = _chromosome_rank(annotation["chrom"])
    tmp = annotation.assign(_chrom_rank=chrom_rank)
    for chrom, group in tmp.sort_values(["_chrom_rank", "chrompos"]).groupby("chrom", sort=False):
        max_pos = int(np.ceil(group["chrompos"].max()))
        if max_pos <= 0:
            continue
        positions = np.arange(bin_size, max_pos + bin_size, bin_size, dtype=np.int64)
        rows.append(pd.DataFrame({"chrom": str(chrom), "chrompos": positions}))
    if not rows:
        raise ValueError("cannot build bins from an empty annotation")
    bins = pd.concat(rows, ignore_index=True)
    chrom_rank = _chromosome_rank(bins["chrom"])
    offsets = _chromosome_offsets_from_positions(bins, chrom_rank)
    bins["abspos"] = bins["chrom"].map(offsets).astype(float) + bins["chrompos"].astype(float)
    bins["_chrom_rank"] = chrom_rank
    return bins.sort_values(["_chrom_rank", "chrompos"]).drop(columns=["_chrom_rank"]).reset_index(drop=True)


def filter_copykat_chromosomes(frame: pd.DataFrame, *, genome: str) -> pd.DataFrame:
    """Match CopyKAT's default human autosome/X chromosome output scope."""

    if genome.lower().startswith("hg"):
        rank = _chromosome_rank(frame["chrom"])
        return frame.loc[rank <= 23].reset_index(drop=True)
    return frame.reset_index(drop=True)


def _chromosome_rank(chrom: pd.Series) -> pd.Series:
    cleaned = chrom.astype(str).str.replace("^chr", "", regex=True)

    def key(value: str) -> int:
        upper = value.upper()
        if upper == "X":
            return 23
        if upper == "Y":
            return 24
        if upper in {"M", "MT"}:
            return 25
        try:
            return int(float(upper))
        except ValueError:
            return 1000 + abs(hash(upper)) % 1000

    return cleaned.map(key)


def _chromosome_offsets(annotation: pd.DataFrame, chrom_rank: pd.Series) -> dict[str, float]:
    tmp = annotation.assign(_chrom_rank=chrom_rank)
    return _chromosome_offsets_from_positions(tmp.rename(columns={"chrompos": "_pos"}), chrom_rank, pos_col="_pos")


def _chromosome_offsets_from_positions(
    positions: pd.DataFrame, chrom_rank: pd.Series, *, pos_col: str = "chrompos"
) -> dict[str, float]:
    tmp = positions.assign(_chrom_rank=chrom_rank)
    offsets: dict[str, float] = {}
    offset = 0.0
    for chrom, group in tmp.sort_values("_chrom_rank").groupby("chrom", sort=False):
        offsets[str(chrom)] = offset
        offset += float(group[pos_col].max()) + 1.0
    return offsets
