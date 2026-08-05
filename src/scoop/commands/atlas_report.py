from __future__ import annotations

import argparse

from scoop.atlas import report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate the scVI annotation quality report and scPoli supervision set."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--input")
    parser.add_argument("--output")
    args = parser.parse_args()
    for path in report(config_path=args.config, input_path=args.input, output_path=args.output):
        print(path)
    return 0
