from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from .storage import append_jsonl


def log_decision(
    run_root: str | Path,
    *,
    module: str,
    decision: str,
    reason: str,
    parameters: Mapping[str, Any] | None = None,
    evidence: Mapping[str, Any] | None = None,
    fallback_used: bool = False,
    fallback_reason: str | None = None,
    human_review_required: bool = False,
    review_reason: str | None = None,
) -> None:
    record = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "module": module,
        "decision": decision,
        "reason": reason,
        "parameters": dict(parameters or {}),
        "evidence": dict(evidence or {}),
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "human_review_required": human_review_required,
        "review_reason": review_reason,
    }
    append_jsonl(record, Path(run_root) / "logs" / "decision_log.jsonl")
