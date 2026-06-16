from __future__ import annotations

import scoop


def test_scoop_package_exports_official_name() -> None:
    assert scoop.OFFICIAL_NAME == "SCOOP: Single Cell Omics Operating Protocol"
