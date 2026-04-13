"""
ARIA Workflow Orchestrator

Coordinates the analysis pipeline, invoking modules based on
decisions from the reasoning engine.
"""

import json
import logging
from pathlib import Path
from typing import Optional

from .reasoning import ReasoningEngine, AnalysisState, DEResult

logger = logging.getLogger(__name__)


class WorkflowOrchestrator:
    """
    Orchestrates the full RNA-seq analysis workflow.

    Unlike traditional fixed pipelines, the orchestrator adapts its
    execution plan based on intermediate results.
    """

    def __init__(self, config: dict, reasoning_engine: Optional[ReasoningEngine] = None):
        self.config = config
        self.engine = reasoning_engine or ReasoningEngine()
        self.execution_log: list[dict] = []
        self.output_dir = Path(config.get("output_dir", "output"))
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self, samplesheet: str):
        """
        Execute the full adaptive analysis pipeline.

        The pipeline adapts at each decision point based on
        intermediate results.
        """
        logger.info("ARIA analysis started")
        self._log_step("pipeline_start", {"samplesheet": samplesheet})

        # Step 1: QC Assessment
        qc_results = self._run_qc(samplesheet)
        qc_decision = self.engine.evaluate_qc(qc_results)
        self._log_step("qc_decision", qc_decision.__dict__)

        if "FAIL" in qc_decision.action:
            logger.error("QC failed. Aborting.")
            return

        # Step 2: Alignment & Quantification
        quant_results = self._run_quantification(samplesheet)
        self._log_step("quantification_complete", {"status": "success"})

        # Step 3: Differential Expression
        de_results = self._run_de_analysis()
        de_decisions = self.engine.evaluate_de_results(de_results)
        for d in de_decisions:
            self._log_step("de_decision", d.__dict__)

        # Step 4: Adaptive pathway analysis
        needs_gsea = any("GSEA" in d.action for d in de_decisions)
        if needs_gsea:
            self._run_gsea()
            self._log_step("gsea_triggered", {"reason": "adaptive_decision"})

        # Step 5: Check for unexpected signatures
        self._check_signatures(de_results)

        # Step 6: Cross-method validation
        self._run_method_comparison()

        # Step 7: Generate report
        self._generate_report()

        logger.info("ARIA analysis complete")
        self._save_execution_log()

    def _run_qc(self, samplesheet):
        """Run QC assessment using MultiQC metrics."""
        logger.info("Running QC assessment...")
        # Implementation: parse MultiQC output or run FastQC/MultiQC
        return []

    def _run_quantification(self, samplesheet):
        """Run STAR + Salmon quantification."""
        logger.info("Running STAR + Salmon quantification...")
        # Implementation: invoke nf-core/rnaseq or custom STAR+Salmon
        return {}

    def _run_de_analysis(self):
        """Run DESeq2 differential expression analysis."""
        logger.info("Running DESeq2 DE analysis...")
        # Implementation: generate and execute R script
        return []

    def _run_gsea(self):
        """Run GSEA with fgsea across all gene set collections."""
        logger.info("Running GSEA (Hallmark, GO:BP/CC/MF, KEGG, Reactome)...")
        # Implementation: generate and execute R script with fgsea

    def _check_signatures(self, de_results):
        """Check for unexpected gene signatures in DE results."""
        logger.info("Checking for unexpected signatures...")
        # Implementation: check for ependymal, immune, etc. markers

    def _run_method_comparison(self):
        """Cross-validate DE results with edgeR and limma-voom."""
        logger.info("Running method comparison...")
        # Implementation: run edgeR + limma-voom

    def _generate_report(self):
        """Generate comprehensive HTML report."""
        logger.info("Generating report...")
        # Implementation: create HTML with embedded figures

    def _log_step(self, step_name: str, details: dict):
        """Log an execution step."""
        entry = {"step": step_name, "details": details}
        self.execution_log.append(entry)
        logger.info(f"Step: {step_name}")

    def _save_execution_log(self):
        """Save the complete execution log."""
        log_path = self.output_dir / "execution_log.json"
        with open(log_path, "w") as f:
            json.dump(self.execution_log, f, indent=2, default=str)
        logger.info(f"Execution log saved to {log_path}")
