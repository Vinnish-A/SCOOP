from __future__ import annotations

import argparse

from scoop.atlas import annotate


def main() -> int:
    parser = argparse.ArgumentParser(description="Produce and apply scVI-space curated atlas labels.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument(
        "--labels", help="curated cluster label TSV/CSV/Parquet; omitted to emit an unresolved template"
    )
    args = parser.parse_args()
    for path in annotate(
        config_path=args.config, input_path=args.input, output_path=args.output, labels_path=args.labels
    ):
        print(path)
    return 0
