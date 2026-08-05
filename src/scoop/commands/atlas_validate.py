from __future__ import annotations

import argparse

from scoop.atlas import validate_atlas


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate atlas embeddings and label-transfer outputs.")
    parser.add_argument("--config", help="run config used for the stage receipt")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(validate_atlas(input_path=args.input, output_path=args.output))
    return 0
