# ARIA: Adaptive Reasoning for Integrated Analysis — An LLM-Powered Framework for Autonomous Transcriptome Analysis with Decision-Aware Workflow Orchestration

## Abstract

RNA-seq transcriptome analysis requires a multi-step workflow involving quality control, alignment, quantification, differential expression testing, pathway analysis, and biological interpretation. While automated pipelines such as nf-core/rnaseq execute these steps reproducibly, the critical decisions between steps — evaluating quality metrics, selecting statistical methods, adapting analysis strategies based on intermediate results, and interpreting findings in biological context — remain dependent on expert bioinformaticians. Here we present ARIA (Adaptive Reasoning for Integrated Analysis), an open-source framework that uses a Large Language Model (LLM) as a reasoning engine to autonomously navigate the decision space of RNA-seq analysis. ARIA implements eight Decision Points (DPs) that govern quality assessment, strategy adaptation, method selection, and result interpretation, combining rule-based thresholds with LLM-driven contextual reasoning. We benchmark ARIA on four public RNA-seq datasets spanning easy (SEQC, ~8,000 DEGs), moderate (Bottomly mouse brain, ~300 DEGs), complex (Airway paired design, ~2,000 DEGs), and hard (Fmr1 KO, ~100 DEGs) scenarios. ARIA correctly identifies paired experimental designs (increasing DEG detection by 21-47%), adaptively switches from DEG-centric to pathway-centric analysis when differential expression is limited, cross-validates results across multiple statistical methods, and generates publication-ready reports with literature-contextualized interpretations. All 7/7 known dexamethasone target genes were recovered in the Airway benchmark. ARIA is freely available at https://github.com/shoo99/ARIA.

## 1. Introduction

### 1.1 The Decision Problem in Transcriptome Analysis

Bulk RNA-seq has become the standard method for profiling transcriptome-wide gene expression changes across biological conditions (Stark et al., 2019). A typical analysis involves multiple sequential steps: quality assessment of raw sequencing data, read alignment to a reference genome, transcript quantification, normalization, differential expression (DE) testing, functional enrichment analysis, and biological interpretation (Conesa et al., 2016). While each step employs well-established tools — STAR for alignment (Dobin et al., 2013), Salmon for quantification (Patro et al., 2017), DESeq2 for DE testing (Love et al., 2014) — the decisions connecting these steps are non-trivial and significantly impact analysis outcomes.

For instance, the number of differentially expressed genes (DEGs) detected in a DE analysis directly determines the downstream analytical strategy: datasets with hundreds of DEGs support traditional over-representation analysis (ORA), while those with few DEGs require Gene Set Enrichment Analysis (GSEA; Subramanian et al., 2005), which leverages the full ranked gene list rather than an arbitrary significance cutoff. Similarly, the presence of batch effects, paired experimental designs, or unexpected cell-type signatures each demands specific analytical responses that are typically recognized and addressed only by experienced bioinformaticians.

### 1.2 Limitations of Existing Automation

Current workflow managers and standardized pipelines address the reproducibility and scalability challenges of RNA-seq analysis. Tools such as Nextflow (Di Tommaso et al., 2017) and Snakemake (Mölder et al., 2021) ensure that predefined analytical steps execute reliably, while community-curated pipelines like nf-core/rnaseq (Ewels et al., 2020; Harshil et al., 2024) provide best-practice default parameters. Web-based platforms such as Galaxy (Afgan et al., 2018) and iDEP (Ge et al., 2018) offer graphical interfaces that lower the barrier to entry.

However, these tools share a fundamental limitation: they execute fixed workflows. The analytical path is determined before data is examined, and adaptation to unexpected results requires manual intervention. When a dataset yields unexpectedly few DEGs, a fixed pipeline continues with ORA — potentially missing biologically meaningful pathway-level changes that GSEA would detect. When a dataset contains an unexpected cell-type signature (e.g., ependymal cell contamination in brain tissue dissections), no automated system flags this or adds compensatory analyses. The expert knowledge required to make these adaptive decisions remains the primary bottleneck in transcriptomic analysis, particularly for laboratories without dedicated bioinformatics support.

### 1.3 Large Language Models as Reasoning Engines

Recent advances in Large Language Models (LLMs) have demonstrated remarkable capabilities in code generation, scientific reasoning, and domain-specific knowledge synthesis (Bubeck et al., 2023). In bioinformatics, LLMs have been applied to literature mining (Luo et al., 2022), variant interpretation (Toufiq et al., 2023), and code assistance for analysis scripts. However, these applications treat the LLM as a passive tool — generating code snippets or answering queries — rather than as an autonomous reasoning agent capable of navigating a complete analytical workflow.

The concept of LLM-based agents that plan, execute, and adapt has gained traction in software engineering (Jimenez et al., 2024) and scientific discovery (Boiko et al., 2023). Yet, to our knowledge, no framework has systematically applied LLM-based adaptive reasoning to the full RNA-seq analysis pipeline, where the agent must not only execute analyses but also evaluate results, recognize patterns, and make informed decisions about subsequent steps.

### 1.4 Our Contribution

We present ARIA (Adaptive Reasoning for Integrated Analysis), an open-source framework that uses an LLM as a reasoning engine to orchestrate RNA-seq analysis with decision-aware workflow management. ARIA's key contributions are:

1. **Decision Protocol with 8 Decision Points (DPs):** A formalized set of decision points covering QC assessment, DE strategy selection, experimental design recognition, unexpected signature detection, cross-method validation, literature-based interpretation, sensitivity analysis, and report generation.

2. **Adaptive Strategy Switching:** Unlike fixed pipelines, ARIA evaluates intermediate results and adapts its analytical strategy. When DEGs are insufficient, it automatically augments the analysis with GSEA, relaxed thresholds, and multi-factor models. When unexpected patterns emerge (e.g., ependymal markers in neuronal tissue), it adds cell-type deconvolution.

3. **Cross-method Validation:** ARIA routinely compares results across DESeq2, edgeR, and limma-voom, providing confidence through concordance rather than relying on a single method.

4. **Automated Contextualization:** The LLM component provides literature-grounded biological interpretation, connecting observed changes to known biology and generating testable hypotheses.

5. **Systematic Benchmarking:** We evaluate ARIA on four public datasets spanning different difficulty levels, demonstrating its adaptive capabilities and quantifying its agreement with expert analyses.

ARIA is designed not to replace bioinformaticians but to augment them — handling routine decisions autonomously while flagging complex situations for expert review. The framework, Docker images, benchmark scripts, and full execution logs are freely available at https://github.com/shoo99/ARIA under the MIT license.

## References

- Afgan, E., et al. (2018). The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update. Nucleic Acids Research, 46(W1), W537-W544.
- Boiko, D. A., et al. (2023). Autonomous chemical research with large language models. Nature, 624, 570-578.
- Bubeck, S., et al. (2023). Sparks of artificial general intelligence: Early experiments with GPT-4. arXiv:2303.12712.
- Conesa, A., et al. (2016). A survey of best practices for RNA-seq data analysis. Genome Biology, 17, 13.
- Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35, 316-319.
- Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15-21.
- Ewels, P. A., et al. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38, 276-278.
- Ge, S. X., et al. (2018). iDEP: an integrated web application for differential expression and pathway analysis of RNA-Seq data. BMC Bioinformatics, 19, 534.
- Jimenez, C. E., et al. (2024). SWE-bench: Can language models resolve real-world GitHub issues? ICLR 2024.
- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
- Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33.
- Patro, R., et al. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14, 417-419.
- Stark, R., Grzelak, M., & Hadfield, J. (2019). RNA sequencing: the teenage years. Nature Reviews Genetics, 20, 631-656.
- Subramanian, A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 102, 15545-15550.
