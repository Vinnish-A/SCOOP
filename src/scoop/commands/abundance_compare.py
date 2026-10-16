"""Compare donor compositions and fit cross-validated binary bag models."""
import argparse

from fastde.abundance_compare import run_comparison
from scsp_agent_sop.config import read_yaml, resolve_run_root


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = read_yaml(args.config)
    root = resolve_run_root(args.config, cfg)
    run_comparison(cfg["abundance"], root / "differential-expression")
