"""Copies of a packaged default must be byte-identical to the one `scoop init` and the commands read."""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
COPIES = {
    "default-run.yaml": ("configs/default_run.yaml", "skill/en/scoop/assets/config/run.yaml",
                         "skill/scoop/assets/config/run.yaml"),
    "annotation-decision-schema.yaml": ("configs/annotation_decision_schema.yaml",
                                        "skill/en/scoop/assets/schemas/annotation-decision.yaml",
                                        "skill/scoop/assets/schemas/annotation-decision.yaml"),
    "h5ad-schema.yaml": ("configs/h5ad_schema.yaml", "skill/en/scoop/assets/schemas/h5ad.yaml",
                         "skill/scoop/assets/schemas/h5ad.yaml"),
}


@pytest.mark.parametrize("packaged,copy", [(name, copy) for name, copies in COPIES.items() for copy in copies])
def test_copy_matches_packaged_default(packaged, copy):
    source = ROOT / "src" / "scoop" / "data" / packaged
    assert (ROOT / copy).read_bytes() == source.read_bytes(), f"{copy} drifted from {source.relative_to(ROOT)}; copy it over"
