from __future__ import annotations

import argparse
import json

from .runtime import detect_capabilities


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="fastcore")
    sub = parser.add_subparsers(dest="cmd")
    sub.add_parser("capabilities")
    args = parser.parse_args(argv)
    if args.cmd in {None, "capabilities"}:
        print(json.dumps(detect_capabilities().to_dict(), indent=2, sort_keys=True))
        return
    raise SystemExit(f"unknown fastcore command: {args.cmd}")
