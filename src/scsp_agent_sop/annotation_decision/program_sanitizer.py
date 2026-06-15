from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class ProgramSanitizerResult:
    program: str
    program_class: str
    action: str
    top_genes: tuple[str, ...]
    matched_genes: tuple[str, ...]
    reason: str

    def to_dict(self) -> dict:
        return asdict(self)


class ProgramSanitizer:
    hemoglobin_prefixes = ("HBA", "HBB", "HBD", "HBE", "HBG", "HBM", "HBQ", "HBZ")
    cell_cycle_genes = frozenset({"MKI67", "TOP2A", "UBE2C", "HMGB2", "PCNA", "TYMS", "CCNB1", "CCNB2", "CDK1", "MCM5"})
    stress_genes = frozenset({"HSPA1A", "HSPA1B", "HSP90AA1", "JUN", "FOS", "DUSP1", "ATF3", "IER2", "DNAJB1"})
    incompatible_lineage_markers = frozenset({"PTPRC", "EPCAM", "KRT8", "KRT18", "PECAM1", "VWF", "COL1A1", "DCN"})

    def __init__(self, *, fraction_threshold: float = 0.35, overlap_threshold: int = 2) -> None:
        self.fraction_threshold = fraction_threshold
        self.overlap_threshold = overlap_threshold

    def classify(self, program: str, top_genes) -> ProgramSanitizerResult:
        genes = tuple(str(gene).upper() for gene in top_genes if str(gene).strip())
        if not genes:
            return ProgramSanitizerResult(program, "unknown", "review", (), (), "no genes supplied")
        checks = [
            ("ribosomal", self._prefix_fraction(genes, ("RPL", "RPS")), "exclude_from_identity", "ribosomal gene prefix enrichment"),
            ("mitochondrial", self._prefix_fraction(genes, ("MT-",)), "exclude_from_identity", "mitochondrial gene prefix enrichment"),
            ("hemoglobin", self._prefix_fraction(genes, self.hemoglobin_prefixes), "exclude_from_identity", "hemoglobin gene enrichment"),
        ]
        for cls, matched, action, reason in checks:
            if len(matched) / len(genes) >= self.fraction_threshold:
                return ProgramSanitizerResult(program, cls, action, genes, matched, reason)
        cell_cycle = tuple(gene for gene in genes if gene in self.cell_cycle_genes)
        if len(cell_cycle) >= self.overlap_threshold:
            return ProgramSanitizerResult(program, "cell_cycle", "state_only", genes, cell_cycle, "cell-cycle marker overlap")
        stress = tuple(gene for gene in genes if gene in self.stress_genes)
        if len(stress) >= self.overlap_threshold:
            return ProgramSanitizerResult(program, "stress", "state_only", genes, stress, "stress/immediate-early marker overlap")
        lineage_hits = tuple(gene for gene in genes if gene in self.incompatible_lineage_markers)
        if len(lineage_hits) >= 3:
            return ProgramSanitizerResult(program, "mixed", "review", genes, lineage_hits, "multiple incompatible lineage markers")
        return ProgramSanitizerResult(program, "identity", "keep_for_identity", genes, (), "no nuisance signature detected")

    def _prefix_fraction(self, genes: tuple[str, ...], prefixes: tuple[str, ...]) -> tuple[str, ...]:
        return tuple(gene for gene in genes if gene.startswith(prefixes))
