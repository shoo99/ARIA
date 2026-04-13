"""Tests for the ARIA reasoning engine."""

import pytest
from aria.core.reasoning import ReasoningEngine, QCMetrics, DEResult


class TestQCDecision:
    """Test Decision Point 1: QC evaluation."""

    def test_all_pass(self):
        engine = ReasoningEngine()
        metrics = [
            QCMetrics("sample1", mapping_rate=0.92, rrna_rate=0.01, five_three_bias=1.1),
            QCMetrics("sample2", mapping_rate=0.90, rrna_rate=0.02, five_three_bias=1.0),
        ]
        decision = engine.evaluate_qc(metrics)
        assert "PASS" in decision.action

    def test_low_mapping(self):
        engine = ReasoningEngine()
        metrics = [
            QCMetrics("sample1", mapping_rate=0.70, rrna_rate=0.01, five_three_bias=1.0),
        ]
        decision = engine.evaluate_qc(metrics)
        assert "WARNING" in decision.action
        assert "low mapping" in decision.rationale

    def test_high_rrna(self):
        engine = ReasoningEngine()
        metrics = [
            QCMetrics("sample1", mapping_rate=0.90, rrna_rate=0.10, five_three_bias=1.0),
        ]
        decision = engine.evaluate_qc(metrics)
        assert "WARNING" in decision.action
        assert "rRNA" in decision.rationale


class TestDEDecision:
    """Test Decision Point 2: DE result evaluation."""

    def test_sufficient_degs(self):
        engine = ReasoningEngine()
        results = [DEResult(tissue="dStr", method="DESeq2", n_deg_lfc1=200)]
        decisions = engine.evaluate_de_results(results)
        assert len(decisions) == 0  # No adaptive action needed

    def test_few_degs_triggers_gsea(self):
        engine = ReasoningEngine()
        results = [DEResult(tissue="Nac", method="DESeq2", n_deg_lfc1=6)]
        decisions = engine.evaluate_de_results(results)
        assert any("GSEA" in d.action for d in decisions)

    def test_very_few_degs_prioritizes_gsea(self):
        engine = ReasoningEngine()
        results = [DEResult(tissue="OFC", method="DESeq2", n_deg_lfc1=3)]
        decisions = engine.evaluate_de_results(results)
        assert any("PRIORITIZE GSEA" in d.action for d in decisions)


class TestDecisionLog:
    """Test decision logging."""

    def test_log_accumulates(self):
        engine = ReasoningEngine()
        engine.evaluate_qc([QCMetrics("s1", mapping_rate=0.9, rrna_rate=0.01, five_three_bias=1.0)])
        engine.evaluate_de_results([DEResult(tissue="dStr", method="DESeq2", n_deg_lfc1=5)])
        log = engine.get_decision_log()
        assert len(log) >= 2
        assert log[0]["decision_point"] == "DP1"
