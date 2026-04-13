"""
ARIA Reasoning Engine

The core LLM-based reasoning engine that evaluates analysis results
and makes decisions about the next analysis steps.
"""

import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class AnalysisState(Enum):
    """Current state of the analysis pipeline."""
    INITIALIZED = "initialized"
    QC_COMPLETE = "qc_complete"
    ALIGNMENT_COMPLETE = "alignment_complete"
    QUANTIFICATION_COMPLETE = "quantification_complete"
    DE_COMPLETE = "de_complete"
    PATHWAY_COMPLETE = "pathway_complete"
    INTERPRETATION_COMPLETE = "interpretation_complete"
    REPORT_GENERATED = "report_generated"


@dataclass
class QCMetrics:
    """Quality control metrics for a sample."""
    sample_name: str
    total_reads: int = 0
    mapping_rate: float = 0.0
    uniquely_mapped_rate: float = 0.0
    rrna_rate: float = 0.0
    duplication_rate: float = 0.0
    five_three_bias: float = 1.0
    gc_content: float = 0.0


@dataclass
class DEResult:
    """Differential expression result summary."""
    tissue: str
    method: str
    n_tested: int = 0
    n_deg_lfc1: int = 0
    n_deg_lfc05: int = 0
    n_deg_nolfc: int = 0
    n_up: int = 0
    n_down: int = 0
    top_genes: list = field(default_factory=list)


@dataclass
class DecisionRecord:
    """Record of a decision made by the reasoning engine."""
    decision_point: str
    trigger: str
    action: str
    rationale: str
    confidence: float = 0.0


class ReasoningEngine:
    """
    LLM-based reasoning engine for adaptive analysis decisions.

    This engine evaluates intermediate results and decides on the next
    analysis steps, mimicking expert bioinformatician decision-making.
    """

    def __init__(self, llm_backend: str = "claude", api_key: Optional[str] = None):
        self.llm_backend = llm_backend
        self.api_key = api_key or os.environ.get("ARIA_API_KEY")
        self.state = AnalysisState.INITIALIZED
        self.decisions: list[DecisionRecord] = []
        self.qc_metrics: list[QCMetrics] = []
        self.de_results: list[DEResult] = []

    def evaluate_qc(self, metrics: list[QCMetrics]) -> DecisionRecord:
        """
        Decision Point 1: Evaluate QC metrics and decide if data is usable.

        Thresholds based on community consensus:
        - Mapping rate >= 85%
        - rRNA < 5%
        - 5'-3' bias: 0.8-1.2
        """
        self.qc_metrics = metrics
        issues = []

        for m in metrics:
            if m.mapping_rate < 0.85:
                issues.append(f"{m.sample_name}: low mapping ({m.mapping_rate:.1%})")
            if m.rrna_rate > 0.05:
                issues.append(f"{m.sample_name}: high rRNA ({m.rrna_rate:.1%})")
            if m.five_three_bias < 0.8 or m.five_three_bias > 1.2:
                issues.append(f"{m.sample_name}: 5'-3' bias ({m.five_three_bias:.2f})")

        if not issues:
            decision = DecisionRecord(
                decision_point="DP1",
                trigger="QC metrics evaluated",
                action="PASS - proceed to DE analysis",
                rationale="All samples pass QC thresholds",
                confidence=0.95
            )
        else:
            decision = DecisionRecord(
                decision_point="DP1",
                trigger="QC metrics evaluated",
                action=f"WARNING - {len(issues)} issues detected",
                rationale="; ".join(issues),
                confidence=0.7
            )

        self.decisions.append(decision)
        self.state = AnalysisState.QC_COMPLETE
        return decision

    def evaluate_de_results(self, results: list[DEResult]) -> list[DecisionRecord]:
        """
        Decision Point 2: Evaluate DE results and decide next strategy.

        If DEGs are insufficient, triggers GSEA and relaxed cutoff analysis.
        """
        self.de_results = results
        decisions = []

        for res in results:
            if res.n_deg_lfc1 < 50:
                # Few DEGs - switch to pathway-level analysis
                decision = DecisionRecord(
                    decision_point="DP2",
                    trigger=f"{res.tissue}: only {res.n_deg_lfc1} DEGs at |LFC|>1",
                    action="ADD GSEA analysis; TRY relaxed LFC (0.5); CONSIDER multi-factor model",
                    rationale="Insufficient DEGs for ORA-based pathway analysis. "
                              "GSEA uses full gene ranking and is more sensitive with few DEGs "
                              "(Subramanian et al., 2005, PNAS).",
                    confidence=0.9
                )
                decisions.append(decision)

            if res.n_deg_lfc1 < 10:
                decision = DecisionRecord(
                    decision_point="DP2",
                    trigger=f"{res.tissue}: very few DEGs ({res.n_deg_lfc1})",
                    action="PRIORITIZE GSEA over ORA; REPORT pathway-level findings as primary result",
                    rationale="With <10 DEGs, individual gene-level analysis is not informative. "
                              "Pathway-level analysis should be the primary analytical framework.",
                    confidence=0.95
                )
                decisions.append(decision)

        self.decisions.extend(decisions)
        self.state = AnalysisState.DE_COMPLETE
        return decisions

    def evaluate_unexpected_signature(self, signature_type: str,
                                       details: str) -> DecisionRecord:
        """
        Decision Point 4: Handle unexpected signatures in the data.
        """
        decision = DecisionRecord(
            decision_point="DP4",
            trigger=f"Unexpected signature detected: {signature_type}",
            action="ADD cell type deconvolution analysis",
            rationale=f"Unexpected {signature_type} signature ({details}) may indicate "
                      "cell type composition differences. Computational deconvolution "
                      "can help distinguish biological signal from technical artifact.",
            confidence=0.8
        )
        self.decisions.append(decision)
        return decision

    def decide_method_comparison(self, primary_method: str,
                                  n_deg_primary: int) -> DecisionRecord:
        """
        Decision Point 5: Decide whether to cross-validate with other methods.
        """
        decision = DecisionRecord(
            decision_point="DP5",
            trigger=f"Primary method ({primary_method}) found {n_deg_primary} DEGs",
            action="RUN edgeR exact test + limma-voom for cross-validation",
            rationale="Cross-method validation increases confidence in DEG calls. "
                      "Genes identified by multiple methods are more robust.",
            confidence=0.85
        )
        self.decisions.append(decision)
        return decision

    def get_decision_log(self) -> list[dict]:
        """Return all decisions as a list of dictionaries."""
        return [
            {
                "decision_point": d.decision_point,
                "trigger": d.trigger,
                "action": d.action,
                "rationale": d.rationale,
                "confidence": d.confidence
            }
            for d in self.decisions
        ]
