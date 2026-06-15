from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def compare_predictions(reference: Path, candidate: Path) -> dict:
    ref = pd.read_csv(reference, sep="\t")
    cand = pd.read_csv(candidate, sep="\t")
    merged = ref.merge(cand, on="cell.names", suffixes=("_reference", "_candidate"))
    if merged.empty:
        raise ValueError("no overlapping cells")
    ref_col = "copykat.pred_reference"
    cand_col = "copykat.pred_candidate"
    confusion = pd.crosstab(merged[ref_col], merged[cand_col])
    return {
        "reference": str(reference),
        "candidate": str(candidate),
        "reference_cells": int(ref.shape[0]),
        "candidate_cells": int(cand.shape[0]),
        "overlap_cells": int(merged.shape[0]),
        "agreement": float((merged[ref_col] == merged[cand_col]).mean()),
        "reference_counts": ref["copykat.pred"].value_counts().to_dict(),
        "candidate_counts": cand["copykat.pred"].value_counts().to_dict(),
        "confusion": {
            str(row): {str(col): int(value) for col, value in confusion.loc[row].items()}
            for row in confusion.index
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare CopyKAT-compatible prediction files.")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    result = compare_predictions(Path(args.reference), Path(args.candidate))
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output_json:
        Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output_json).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

