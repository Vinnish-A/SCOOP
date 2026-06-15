from __future__ import annotations

from scsp_agent_sop.annotation_decision.program_sanitizer import ProgramSanitizer


def test_ribosomal_program_classifies_correctly() -> None:
    result = ProgramSanitizer().classify("P1", ["RPL3", "RPS6", "RPLP0", "ACTB"])
    assert result.program_class == "ribosomal"
    assert result.action == "exclude_from_identity"


def test_cell_cycle_program_classifies_correctly() -> None:
    result = ProgramSanitizer().classify("P2", ["MKI67", "TOP2A", "UBE2C", "PCNA"])
    assert result.program_class == "cell_cycle"
    assert result.action == "state_only"


def test_unknown_program_returns_identity_or_review() -> None:
    result = ProgramSanitizer().classify("P3", ["GENE1", "GENE2"])
    assert result.program_class in {"identity", "unknown", "mixed"}
