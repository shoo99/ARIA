---
title: 'ARIA: Adaptive Reasoning for Integrated Analysis of RNA-seq Data'
tags:
  - Python
  - bioinformatics
  - RNA-seq
  - transcriptomics
  - large language models
  - differential expression
authors:
  - name: Byeongsoo Kang
    orcid: 0009-0007-2324-2351
    affiliation: 1
affiliations:
  - name: SYSOFT, Republic of Korea
    index: 1
date: 21 April 2026
bibliography: paper.bib
---

# Summary

ARIA (Adaptive Reasoning for Integrated Analysis) is a Python framework that uses a Large Language Model (LLM) as a reasoning engine to autonomously navigate the decision space of RNA-seq transcriptome analysis. Unlike fixed pipelines that execute predetermined workflows regardless of data characteristics, ARIA implements eight Decision Points (DPs) that evaluate intermediate results and adaptively select downstream analysis strategies — mimicking how an experienced bioinformatician makes context-dependent decisions.

# Statement of Need

RNA-seq analysis requires numerous decisions between computational steps: evaluating quality metrics, selecting statistical methods, recognizing experimental designs, and interpreting results in biological context. While workflow managers like Nextflow and platforms like Galaxy automate execution, the critical reasoning that connects steps remains dependent on expert bioinformaticians. This creates a bottleneck in laboratories without dedicated computational staff and introduces variability when different analysts make different choices on the same data.

ARIA addresses this gap by formalizing the decision-making process into a structured protocol. For example, when a dataset yields fewer than 100 differentially expressed genes, ARIA automatically switches from over-representation analysis to Gene Set Enrichment Analysis (GSEA), which leverages the full ranked gene list. When paired experimental designs are detected in metadata, ARIA adjusts the statistical model to account for within-subject correlation, increasing statistical power by 21–47% in validated benchmarks.

# Key Features

- **Eight Decision Points** covering QC evaluation (DP1), DE strategy selection (DP2), experimental design recognition (DP3), signature detection (DP4), cross-method validation (DP5), literature interpretation (DP6), sensitivity analysis (DP7), and automated reporting (DP8).
- **Adaptive strategy switching** based on intermediate results: the analysis path is determined by the data, not by a fixed configuration.
- **Multi-method cross-validation** using DESeq2, edgeR, and limma-voom with concordance assessment (r > 0.99 across benchmarks).
- **Docker-based execution** ensuring full computational reproducibility of R-based statistical analyses.
- **Systematic benchmarking** on four public datasets spanning three species (human, mouse, *Drosophila*) with varying complexity.

# Benchmarking

ARIA has been validated on four public RNA-seq datasets: SEQC (GSE49712, 10,430 DEGs), Airway (GSE52778, 951 DEGs with paired design), Fmr1 KO (GSE180135, 398 DEGs), and Pasilla (GSE18508, 224 DEGs with mixed library types). Key validated behaviors include paired design detection (+21–47% DEG gain), technical covariate adjustment (+4–30%), and recovery of known biological markers (15/16 tissue markers in SEQC, 7/7 dexamethasone targets in Airway, 9/11 FMRP targets in Fmr1 KO).

# Related Preprint

A detailed methodological description, including full benchmark results, ablation studies, reproducibility tests, and metadata obfuscation experiments, is available as a preprint:

> Kang B. (2026). *ARIA: Adaptive Reasoning for Integrated Analysis — An LLM-Powered Framework for Autonomous Transcriptome Analysis with Decision-Aware Workflow Orchestration.* Research Square. [doi:10.21203/rs.3.rs-9500973/v1](https://doi.org/10.21203/rs.3.rs-9500973/v1)

# Acknowledgements

This work was conducted independently at SYSOFT. The author thanks the developers of DESeq2, edgeR, limma, fgsea, and the Anthropic API for the tools that underpin ARIA's analysis capabilities.

# References
