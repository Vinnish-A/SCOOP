from __future__ import annotations

from fastcore.backends.omicverse_cpu import run_omicverse_cpu_core as legacy_cpu
from fastcore.backends.omicverse_mixed import run_omicverse_mixed_core as legacy_mixed
from fastcore.backends.omicverse_rust_oom import run_omicverse_rust_oom_core as legacy_oom
from omicverse_transfer import (
    run_omicverse_cpu_core,
    run_omicverse_mixed_core,
    run_omicverse_rust_oom_core,
)
from omicverse_transfer.vendor.omicverse_gpl import run_vendored_omicverse_cpu_core
from fastcore.vendor.omicverse_gpl import run_vendored_omicverse_cpu_core as legacy_vendored_cpu


def test_omicverse_transfer_exports_core_backends() -> None:
    assert run_omicverse_cpu_core is legacy_cpu
    assert run_omicverse_mixed_core is legacy_mixed
    assert run_omicverse_rust_oom_core is legacy_oom


def test_omicverse_transfer_keeps_vendored_compat_path() -> None:
    assert run_vendored_omicverse_cpu_core is legacy_vendored_cpu
