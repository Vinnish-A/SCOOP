from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .io import write_json


@dataclass
class AbundanceResult:
    mode: str
    results: pd.DataFrame
    predictions: pd.DataFrame
    metrics: dict[str, Any]
    manifest: dict[str, Any]

    def write(self, output_dir: str | Path) -> dict[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        results_path = out / f"abundance_{self.mode}_results.tsv"
        predictions_path = out / f"abundance_{self.mode}_predictions.tsv"
        metrics_path = out / f"abundance_{self.mode}_metrics.json"
        manifest_path = out / "abundance_manifest.json"
        self.results.to_csv(results_path, sep="\t", index=False)
        self.predictions.to_csv(predictions_path, sep="\t", index=False)
        write_json(self.metrics, metrics_path)
        self.manifest.setdefault("outputs", {})
        self.manifest["outputs"].update(
            {
                "results": str(results_path),
                "predictions": str(predictions_path),
                "metrics": str(metrics_path),
                "manifest": str(manifest_path),
            }
        )
        write_json(self.manifest, manifest_path)
        return {
            "results": str(results_path),
            "predictions": str(predictions_path),
            "metrics": str(metrics_path),
            "manifest": str(manifest_path),
        }
