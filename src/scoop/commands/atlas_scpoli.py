from __future__ import annotations

import argparse

from scoop.atlas import train_scpoli


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Train scPoli only after verifying the hash-bound annotation approval."
    )
    parser.add_argument("--config", help="run config used for the stage receipt")
    parser.add_argument("--gate", required=True)
    parser.add_argument("--output")
    parser.add_argument("--epochs", type=int, default=100)
    args = parser.parse_args()
    print(train_scpoli(gate_path=args.gate, output_path=args.output, n_epochs=args.epochs))
    return 0
