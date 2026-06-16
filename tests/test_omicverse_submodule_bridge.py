from __future__ import annotations

import sys
import types

from fastcore.backend_plan import canonical_backend
from fastcore.backends.omicverse_cpu import run_omicverse_cpu_core
from fastcore.vendor.omicverse_gpl import run_vendored_omicverse_cpu_core
from scsp_agent_sop.omicverse_facilities import require_omicverse


def test_legacy_backend_names_are_aliases() -> None:
    assert canonical_backend("omicverse_cpu") == "fastcore_cpu"
    assert canonical_backend("omicverse_cpu_gpu_mixed") == "fastcore_mixed"
    assert canonical_backend("omicverse_rust_oom") == "fastcore_oom"


def test_fastcore_cpu_uses_vendored_compat_path() -> None:
    assert run_omicverse_cpu_core.__name__ == "run_omicverse_cpu_core"
    assert callable(run_vendored_omicverse_cpu_core)


def test_require_omicverse_prefers_imported_module(monkeypatch) -> None:
    fake = types.SimpleNamespace(marker="installed_or_submodule")
    monkeypatch.setitem(sys.modules, "omicverse", fake)
    assert require_omicverse() is fake
