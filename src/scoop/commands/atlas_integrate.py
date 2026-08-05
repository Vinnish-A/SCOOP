from __future__ import annotations

import argparse

from scoop.atlas import integrate


def main() -> int:
    parser = argparse.ArgumentParser(description="Run scVI integration and clustering in X_scVI.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--max-epochs", type=int)
    args = parser.parse_args()
    print(
        integrate(config_path=args.config, input_path=args.input, output_path=args.output, max_epochs=args.max_epochs)
    )
    return 0
