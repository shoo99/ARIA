# ARIA: Adaptive Reasoning for Integrated Analysis

**An LLM-Powered Framework for Autonomous Transcriptome Analysis with Decision-Aware Workflow Orchestration**

**Byeongsoo Kang**¹*

¹ SYSOFT, Republic of Korea

\* Corresponding author: shoo99@gmail.com

---

## Abstract

RNA-seq transcriptome analysis requires a multi-step workflow involving quality control, alignment, quantification, differential expression testing, pathway analysis, and biological interpretation. While automated pipelines such as nf-core/rnaseq execute these steps reproducibly, the critical decisions between steps — evaluating quality metrics, selecting statistical methods, adapting analysis strategies based on intermediate results, and interpreting findings in biological context — remain dependent on expert bioinformaticians. Here we present ARIA (Adaptive Reasoning for Integrated Analysis), an open-source framework that uses a Large Language Model (LLM) as a reasoning engine to autonomously navigate the decision space of RNA-seq analysis. ARIA implements eight Decision Points (DPs) that govern quality assessment, strategy adaptation, method selection, and result interpretation, combining rule-based thresholds with LLM-driven contextual reasoning. We benchmark ARIA on four public RNA-seq datasets spanning three species: SEQC (GSE49712, human, 10,430 DEGs), Airway (GSE52778, human paired design, 951 DEGs), Fmr1 KO (GSE180135, mouse, 398 DEGs), and Pasilla (GSE18508, Drosophila, 224 DEGs with mixed library types). ARIA correctly identifies paired experimental designs (increasing DEG detection by 21–47%), detects technical covariates such as library type (+4–30% DEG gain), adaptively selects analysis strategies based on DEG counts, and cross-validates results across DESeq2/edgeR/limma-voom (r > 0.99). All 7/7 known dexamethasone targets were recovered in Airway, 9/11 FMRP targets in Fmr1 KO, and 15/16 brain/cancer markers in SEQC were correctly assigned. ARIA is freely available at https://github.com/shoo99/ARIA.

**Keywords:** RNA-seq, differential expression, LLM, autonomous analysis, decision-aware pipeline, GSEA

---

## 1. Introduction

### 1.1 The Decision Problem in Transcriptome Analysis

Bulk RNA-seq has become the standard method for profiling transcriptome-wide gene expression changes across biological conditions (Stark et al., 2019). A typical analysis involves multiple sequential steps: quality assessment of raw sequencing data, read alignment to a reference genome, transcript quantification, normalization, differential expression (DE) testing, functional enrichment analysis, and biological interpretation (Conesa et al., 2016). While each step employs well-established tools — STAR for alignment (Dobin et al., 2013), Salmon for quantification (Patro et al., 2017), DESeq2 for DE testing (Love et al., 2014) — the decisions connecting these steps are non-trivial and significantly impact analysis outcomes.

For instance, the number of differentially expressed genes (DEGs) detected in a DE analysis directly determines the downstream analytical strategy: datasets with hundreds of DEGs support traditional over-representation analysis (ORA), while those with few DEGs require Gene Set Enrichment Analysis (GSEA; Subramanian et al., 2005), which leverages the full ranked gene list rather than an arbitrary significance cutoff. Similarly, the presence of batch effects, paired experimental designs, or unexpected cell-type signatures each demands specific analytical responses that are typically recognized and addressed only by experienced bioinformaticians.

### 1.2 Limitations of Existing Automation

Current workflow managers and standardized pipelines address the reproducibility and scalability challenges of RNA-seq analysis. Tools such as Nextflow (Di Tommaso et al., 2017) and Snakemake (Mölder et al., 2021) ensure that predefined analytical steps execute reliably, while community-curated pipelines like nf-core/rnaseq (Ewels et al., 2020) provide best-practice default parameters. Web-based platforms such as Galaxy (Afgan et al., 2018) and iDEP (Ge et al., 2018) offer graphical interfaces that lower the barrier to entry.

However, these tools share a fundamental limitation: they execute fixed workflows. The analytical path is determined before data is examined, and adaptation to unexpected results requires manual intervention. When a dataset yields unexpectedly few DEGs, a fixed pipeline continues with ORA — potentially missing biologically meaningful pathway-level changes that GSEA would detect. The expert knowledge required to make these adaptive decisions remains the primary bottleneck in transcriptomic analysis.

### 1.3 Large Language Models as Reasoning Engines

Recent advances in Large Language Models (LLMs) have demonstrated capabilities in code generation, scientific reasoning, and domain-specific knowledge synthesis (Bubeck et al., 2023). In bioinformatics, LLMs have been applied to literature mining and code assistance. However, these applications treat the LLM as a passive tool rather than as an autonomous reasoning agent capable of navigating a complete analytical workflow.

The concept of LLM-based agents that plan, execute, and adapt has gained traction in software engineering (Jimenez et al., 2024) and scientific discovery (Boiko et al., 2023). Yet, to our knowledge, no framework has systematically applied LLM-based adaptive reasoning to the full RNA-seq analysis pipeline.

### 1.4 Our Contribution

We present ARIA, an open-source framework with five key contributions:

1. **Decision Protocol with 8 Decision Points (DPs)** covering QC, DE strategy, design recognition, signature detection, cross-validation, interpretation, sensitivity analysis, and reporting.
2. **Adaptive Strategy Switching** based on intermediate results.
3. **Cross-method Validation** across DESeq2, edgeR, and limma-voom.
4. **Automated Contextualization** via LLM-driven biological interpretation.
5. **Systematic Benchmarking** on four public datasets of varying difficulty.

---

## 2. Methods

### 2.1 Architecture

ARIA is organized into four layers (Figure 1): (1) **Reasoning Engine** — an LLM serving as the reasoning backbone; (2) **Decision Protocol** — eight formalized Decision Points combining rule-based thresholds with LLM reasoning; (3) **Execution Engine** — Docker-based R/Python script generation and execution; (4) **Analysis Modules** — DESeq2, fgsea, STRING DB, WGCNA, cell type deconvolution, and ssGSEA.

### 2.2 Decision Protocol

Each Decision Point (DP) operates through a hybrid mechanism: rule-based thresholds handle deterministic criteria, while the LLM reasoning engine handles contextual evaluation, edge cases, and interpretive tasks. The table below summarizes each DP, followed by detailed descriptions of the LLM's role in the more complex decision points (DP5–DP7).

| DP | Trigger | Action | Validation |
|----|---------|--------|------------|
| DP1 | QC metrics available | Pass/Warn/Fail based on mapping rate ≥85%, rRNA <5%, 5'–3' bias 0.8–1.2 | All 4 benchmarks |
| DP2 | DEG count evaluated | ≥100: standard; 10–100: add GSEA; <10: prioritize GSEA | SEQC (standard), Airway (standard), Fmr1 (standard) |
| DP3 | Metadata examined | Detect paired/blocked designs, select model formula | Airway (+21–47% DEGs), Pasilla (+4–30%) |
| DP4 | Unexpected gene clusters | Flag and add cell type deconvolution | Demonstrated in real-world use |
| DP5 | Primary DE complete | Cross-validate with edgeR + limma-voom | Fmr1 (r=0.9999), Pasilla (r=0.9906) |
| DP6 | Pathways identified | Literature-based hypothesis generation | Demonstrated in real-world use |
| DP7 | Results available | Sensitivity analysis across LFC cutoffs | All benchmarks |
| DP8 | Analysis complete | Generate HTML report + Excel + figures | All benchmarks |

#### 2.2.1 DP5: Cross-Method Validation — LLM Role

At DP5, the Execution Engine runs DESeq2, edgeR (exact test and quasi-likelihood), and limma-voom on identical data. The LLM's role is to *evaluate concordance* and *diagnose discrepancies*:

1. **Input to LLM:** A structured summary containing (a) DEG counts per method at multiple cutoffs, (b) log2FC Pearson correlation, (c) Venn diagram overlap statistics.
2. **LLM reasoning process:** The LLM receives a prompt: *"Given these cross-method results, assess whether the methods are concordant. If DEG counts differ by >50% or LFC correlation is <0.95, identify which method is the outlier and hypothesize why."*
3. **Decision output:** If concordance is high (r > 0.95, overlap > 70%), the LLM confirms the primary DESeq2 results. If discrepancy is detected, the LLM generates a diagnostic report identifying the likely cause (e.g., edgeR QLF being overly conservative at small n, as observed in our real-world testing where QLF detected only 1–2 DEGs with n=3).
4. **Rule component:** The correlation threshold (r > 0.95) and overlap threshold (>70%) are rule-based; the diagnostic interpretation is LLM-driven.

#### 2.2.2 DP6: Literature-Based Interpretation — LLM Role

DP6 represents the most LLM-dependent decision point. The process operates as follows:

1. **Input to LLM:** (a) Top 10 DEGs with gene names, fold changes, and significance values; (b) Top enriched pathways from GSEA; (c) Experimental context (species, tissue, condition).
2. **Knowledge source:** The LLM draws on its pre-trained knowledge base, which encompasses biomedical literature up to its training cutoff. No external database queries are made at this step — the LLM functions as a compressed literature index.
3. **Prompt structure:** *"Given these DEGs [list] and enriched pathways [list] in [tissue] of [condition] mice, identify: (a) known functions of the top DEGs, (b) whether the enriched pathways are consistent with the known biology of this model, (c) testable hypotheses arising from unexpected findings. Label all interpretations as hypothesis-level."*
4. **Output:** The LLM generates structured annotations for each top DEG and pathway, explicitly labeling confidence levels. Interpretations that align with well-established biology (e.g., Fmr1 downregulation in Fmr1 KO) receive higher confidence than novel associations.
5. **Critical safeguard:** All LLM-generated interpretations are tagged with a disclaimer: *"Hypothesis-level interpretation based on LLM knowledge. Verify citations independently before publication."* This addresses the hallucination risk discussed in Section 4.4.

#### 2.2.3 DP7: Sensitivity Analysis — LLM Role

DP7 is primarily rule-based in execution but LLM-assisted in interpretation:

1. **Execution (rule-based):** The system automatically tests a matrix of cutoff combinations: padj thresholds (0.01, 0.05, 0.1) × LFC thresholds (0, 0.5, 1.0, 2.0), generating a 12-cell sensitivity table.
2. **LLM interpretation:** The LLM receives the sensitivity table and is prompted: *"Analyze this sensitivity table. Identify the cutoff combination that balances biological significance with statistical rigor for this dataset. Consider the number of DEGs, the typical fold-change magnitude, and the sample size."*
3. **Recommendation output:** The LLM recommends an optimal cutoff (e.g., *"With n=3 and subtle KO effects, |LFC|>0.5 at padj<0.05 provides 635 DEGs — sufficient for pathway analysis while maintaining biological relevance"*) with explicit reasoning.

### 2.3 Implementation

ARIA is implemented in Python 3.10+ with R analysis modules executed in Docker containers. Key tools: DESeq2 v1.46.0, edgeR, limma, fgsea, msigdbr, WGCNA, STRING DB v12. All R scripts are generated from parameterized templates for reproducibility.

### 2.4 Benchmark Datasets and Difficulty Classification

We selected four datasets spanning a range of analytical challenges. Difficulty was classified based on three objective criteria: (1) **expected effect size** (large inter-group differences vs. subtle changes), (2) **design complexity** (simple two-group vs. paired/blocked with covariates), and (3) **sample size** relative to biological variability.

| Dataset | GEO | Species | Design | N | Difficulty | Rationale |
|---------|-----|---------|--------|---|------------|-----------|
| SEQC | GSE49712 | Human | UHRR vs HBRR | 10 | Easy | Massive effect size (~10,000 DEGs); reference RNAs with minimal biological variability; no covariates |
| Airway | GSE52778 | Human | Dex ± ASM, paired | 8 | Complex | Moderate effect size; paired design requiring covariate modeling; 4 cell line backgrounds |
| Fmr1 KO | GSE180135 | Mouse | KO vs WT neurons | 6 | Moderate | Moderate effect size; simple two-group design but small n=3 per group; in vitro neurons |
| Pasilla | GSE18508 | Drosophila | pasilla KD vs control | 7 | Moderate | Moderate effect size; mixed SE/PE library types requiring technical covariate; 4 vs 3 unbalanced design |

The "Easy" classification for SEQC reflects the comparison of fundamentally different RNA populations (cancer cell lines vs. brain tissue), producing thousands of DEGs even with conservative thresholds. "Complex" for Airway reflects not the DEG count but the paired experimental design, which requires correct covariate modeling — a decision that substantially impacts results (Table 2). "Moderate" for both Fmr1 KO and Pasilla reflects moderate effect sizes combined with design features that test specific ARIA capabilities (small sample size and technical covariates, respectively).

Count matrices were obtained from GEO supplementary data (SEQC: HTSeq counts; Fmr1: DESeq2 counts) or Bioconductor packages (Airway: `airway` R package; Pasilla: `pasilla` R package).

### 2.5 Evaluation Metrics

(A) DEG accuracy vs known biology (target gene recovery rate); (B) Decision appropriateness (correct strategy selection); (C) Cross-method concordance (LFC correlation); (D) GSEA pathway yield; (E) Analysis time.

---

## 3. Results

### 3.1 Benchmark Performance Overview

**Table 1a. Benchmark datasets and differential expression results**

| | SEQC | Airway | Fmr1 KO | Pasilla |
|---|---|---|---|---|
| **GEO** | GSE49712 | GSE52778 | GSE180135 | GSE18508 |
| **Species** | Human | Human | Mouse | Drosophila |
| **Design** | UHRR vs HBRR | Dex ± (paired) | KO vs WT | KD vs control |
| **Samples** | 5 vs 5 | 4 vs 4 | 3 vs 3 | 4 vs 3 |
| **Difficulty** | Easy | Complex | Moderate | Moderate |
| **Genes tested** | 16,417 | ~15,000 | ~17,000 | 9,686 |
| **DEGs (\|LFC\|>1)** | **10,430** | **951** | **398** | **224** |
| **DEGs (\|LFC\|>0.5)** | 13,818 | 2,426 | 1,654 | 635 |

**Table 1b. Pathway enrichment and validation results**

| | SEQC | Airway | Fmr1 KO | Pasilla |
|---|---|---|---|---|
| **GSEA Hallmark** | 32 | 15 | 17 | — |
| **GSEA GO:BP** | 1,352 | 161 | 607 | — |
| **GSEA Reactome** | 505 | 17 | 198 | — |
| **Known targets** | Brain 9/9, Cancer 6/7 | Dex **7/7** | FMRP **9/11** | Splicing KD |
| **Key DP validated** | DP1, DP2 | DP1, DP2, **DP3** | DP1, DP2, DP5 | DP1, DP2, DP3, DP5, DP7 |
| **Cross-method r** | — | — | 0.9999 | 0.9906 |

### 3.2 DP3: Experimental Design Recognition (Airway)

ARIA's DP3 correctly identified the paired design and selected `~ cell + dex` over naive `~ dex`.

**Table 2. Impact of design recognition**

| Cutoff | Unpaired | Paired (ARIA) | Gain |
|--------|----------|---------------|------|
| padj<0.05, \|LFC\|>1.0 | 785 | **951** | **+21%** |
| padj<0.05, \|LFC\|>0.5 | 1,860 | **2,426** | **+30%** |
| padj<0.05, no LFC | 2,773 | **4,081** | **+47%** |

All 7 known dexamethasone targets (DUSP1, KLF15, CRISPLD2, PER1, FKBP5, TSC22D3, ZBTB16) were recovered with high significance (padj < 1e-19, LFC 2.6–7.4).

### 3.3 Known Target Validation Across Benchmarks

**SEQC (UHRR vs HBRR):** 9/9 brain-enriched genes (GFAP, MBP, SYN1, GAD1, SLC17A7, NEFL, NEFM, NEFH, ENO2) correctly identified as upregulated in HBRR (LFC +2.2 to +12.7). 6/7 cancer genes (MYC, CCND1, CDK4, EGFR, ERBB2, TP53) correctly identified as upregulated in UHRR. BCL2 showed unexpected direction, consistent with known brain expression.

**Fmr1 KO:** Fmr1 itself detected as the top DEG (LFC = −0.98, padj = 3.7e-83). Key FMRP-regulated synaptic proteins detected: Dlg4/PSD-95 (LFC = −0.36), Shank3 (LFC = −0.45), Camk2a (LFC = −0.42). NMDA receptor subunits: Grin1 (LFC = −0.92), Grin2a (LFC = −0.81), Grin2b (LFC = −0.48). All directions consistent with loss of FMRP-mediated translational regulation.

### 3.4 DP5: Cross-Method Concordance

| Benchmark | DESeq2 | edgeR exact | limma-voom | LFC correlation |
|-----------|--------|-------------|------------|-----------------|
| Fmr1 KO | 1,654 | ~1,600 | ~1,500 | **r = 0.9999** |
| Pasilla | 635 | 739 | 674 | **r = 0.9906** |

### 3.5 DP2: Adaptive Strategy Selection

ARIA correctly classified each dataset's difficulty level:
- **SEQC (10,430 DEGs):** → "standard" (ORA + GSEA) ✓
- **Airway (951 DEGs):** → "standard" ✓
- **Fmr1 KO (398 DEGs):** → "standard" ✓
- **Pasilla (224 DEGs):** → "standard" ✓ (also detected SE/PE library type covariate via DP3)

### 3.6 DP3 + DP7: Covariate Detection and Sensitivity Analysis (Pasilla)

The Pasilla dataset (GSE18508) contains mixed single-end and paired-end libraries, creating a technical covariate. ARIA's DP3 detected this library type difference and included it as a covariate (`~ type + condition`), improving DEG detection:

**Table 4. Impact of library type covariate (Pasilla)**

| Cutoff | No covariate | With type (ARIA) | Gain |
|--------|-------------|-----------------|------|
| padj<0.05, \|LFC\|>1.0 | 216 | **224** | **+4%** |
| padj<0.05, \|LFC\|>0.5 | 583 | **635** | **+9%** |
| padj<0.05, no LFC | 853 | **1,108** | **+30%** |

DP7 (sensitivity analysis) generated a 12-combination cutoff table (3 padj thresholds × 4 LFC thresholds), providing a comprehensive view of result stability. Cross-method validation (DP5) showed high concordance: DESeq2 (635), edgeR (739), limma-voom (674), with LFC correlation r = 0.9906.

### 3.7 GSEA Demonstrates Value Beyond DEG Counting

Even in datasets with abundant DEGs, GSEA provided complementary insights. In Fmr1 KO, 607 GO:BP pathways were enriched (padj < 0.05), revealing systemic disruption of synaptic signaling, translation regulation, and neuronal development — biological processes well-documented in Fragile X syndrome literature.

---

## 4. Discussion

### 4.1 Adaptive Decision Making as Core Innovation

ARIA's central contribution is the reasoning layer that adaptively connects established analysis tools. In the Airway benchmark, design recognition (DP3) increased DEG detection by 21–47%. This represents over a thousand genes that would be missed by a naive analysis — a quantifiable improvement from a single adaptive decision.

### 4.2 Comparison with Existing Tools

| Feature | nf-core | Galaxy | iDEP | ARIA |
|---------|---------|--------|------|------|
| Execution automation | ✓ | ✓ | ✓ | ✓ |
| Result-based adaptation | ✗ | ✗ | △ | **✓** |
| Design recognition | ✗ | ✗ | ✗ | **✓** |
| Cross-method validation | ✗ | ✗ | ✗ | **✓** |
| Biological interpretation | ✗ | ✗ | △ | **✓** |
| Decision transparency | ✗ | ✗ | ✗ | **✓** |
| Reproducibility | ✓ | ✓ | △ | ✓ |

**Decision transparency** is defined here as the ability of the system to record, explain, and expose the reasoning behind each analytical choice made during the workflow. This encompasses three levels:

**(1) What was decided:** A structured log of every decision point triggered during the analysis, including the action taken (e.g., "switched to GSEA-priority strategy") and the threshold or criterion that triggered it (e.g., "DEG count = 47, below threshold of 100").

**(2) Why it was decided:** A natural-language rationale for each decision, explaining the reasoning chain. For example: *"With only 47 DEGs at |LFC|>1 and n=3 per group, ORA-based pathway analysis is statistically underpowered. GSEA, which uses the full ranked gene list, is more appropriate for this dataset (cf. Subramanian et al., 2005)."*

**(3) What alternatives were considered:** Where applicable, the log records alternative actions that were evaluated and rejected, with reasons. For example: *"Multi-factor model considered but not applied — single tissue analysis, no blocking factor detected."*

Existing tools do not provide this level of transparency. **nf-core/rnaseq** logs tool versions and execution parameters but does not record analytical reasoning — it follows a fixed path with no decisions to document. **Galaxy** records the history of tool invocations with parameters, but the *why* behind parameter choices is not captured, as these are made by the user through the GUI. **iDEP** provides interactive visualizations that implicitly guide user decisions, but the decision process itself is not logged or exportable.

ARIA's decision log is saved as a structured JSON file (`execution_log.json`) alongside every analysis, enabling: (a) post-hoc audit by reviewers or collaborators, (b) reproducibility of the analytical reasoning (not just the execution), and (c) identification of points where the LLM's reasoning may need expert correction.

### 4.3 Human-in-the-Loop Design

ARIA is designed to augment, not replace, bioinformaticians. Biological interpretations are labeled as hypothesis-level. The decision log enables audit and review. Complex situations are flagged for expert judgment.

### 4.4 Limitations

#### 4.4.1 LLM Hallucination Risk and Mitigation

LLMs can generate plausible but factually incorrect statements — a phenomenon known as "hallucination." In ARIA, this risk is most acute at DP6 (literature-based interpretation), where the LLM may: (a) cite non-existent papers, (b) attribute incorrect functions to genes, or (c) fabricate mechanistic connections.

We implement four mitigation strategies:

**(1) Scope restriction:** The LLM is only asked to interpret genes and pathways that are statistically significant. It does not generate hypotheses from non-significant results, reducing the space for fabrication.

**(2) Explicit labeling:** All LLM-generated interpretations are tagged with `[Hypothesis — verify independently]`. The decision log records the exact prompt and response, enabling post-hoc audit.

**(3) Cross-referencing with rule-based outputs:** LLM interpretations are constrained by the statistical results. If the LLM claims a gene is upregulated but the DE results show downregulation, this contradiction is automatically flagged.

**(4) Human review checkpoints:** ARIA's report includes a dedicated "Interpretations requiring expert review" section that aggregates all LLM-generated biological claims, making it efficient for domain experts to audit.

**Observed hallucination in practice:** During real-world testing, we observed that the LLM occasionally described gene functions with minor inaccuracies (e.g., attributing a protein to the wrong subcellular compartment) and once generated a reference that could not be verified. These instances were caught during expert review, reinforcing the necessity of the human-in-the-loop design.

#### 4.4.2 Statistical Power

Statistical power is constrained by sample size, not by the analytical method. ARIA transparently communicates this limitation and adapts its strategy (e.g., switching to GSEA), but it cannot overcome fundamental statistical constraints inherent to small-n experiments.

#### 4.4.3 Cell Type Confounding

Bulk RNA-seq cannot distinguish gene expression changes within a cell type from changes in cell type composition. ARIA's marker-based deconvolution (DP4) provides an approximate assessment, but single-cell or single-nucleus RNA-seq is required for definitive cell-type-resolved analysis.

#### 4.4.4 API Cost and Reproducibility

ARIA requires LLM API access (~$2–10 per analysis). LLM outputs are inherently stochastic; while statistical analyses use fixed seeds, the interpretive component (DP6) may vary slightly across runs. The decision log mitigates this by recording the exact reasoning for each run.

#### 4.4.5 Scope

The current version of ARIA is designed for bulk RNA-seq two-group comparisons. Extensions to time-series designs, single-cell RNA-seq, and multi-omics integration are planned but not yet implemented.

### 4.5 Future Directions

Extensions to single-cell RNA-seq, multi-omics integration, time-series designs, and domain-specific fine-tuned LLMs are planned.

---

## 5. Data and Code Availability

- **Source code:** https://github.com/shoo99/ARIA
- **License:** MIT
- **Docker:** aria-bench:latest
- **Benchmarks:** GSE49712, GSE52778, GSE180135, GSE18508

---

## Author Contributions

B.K. conceived the project, designed the framework architecture, implemented the software, performed all benchmark analyses, and wrote the manuscript.

## Acknowledgments

The authors thank the developers of DESeq2, edgeR, fgsea, nf-core, and the broader Bioconductor community for making their tools freely available. The SEQC/MAQC-III Consortium, Himes et al., Brooks et al., and Bhatt et al. are acknowledged for making their RNA-seq datasets publicly accessible.

## Competing Interests

The author declares no competing interests.

## References

1. Afgan, E., et al. (2018). Nucleic Acids Research, 46, W537–W544.
2. Boiko, D.A., et al. (2023). Nature, 624, 570–578.
3. Bubeck, S., et al. (2023). arXiv:2303.12712.
4. Conesa, A., et al. (2016). Genome Biology, 17, 13.
5. Di Tommaso, P., et al. (2017). Nature Biotechnology, 35, 316–319.
6. Dobin, A., et al. (2013). Bioinformatics, 29, 15–21.
7. Ewels, P.A., et al. (2020). Nature Biotechnology, 38, 276–278.
8. Ge, S.X., et al. (2018). BMC Bioinformatics, 19, 534.
9. Jimenez, C.E., et al. (2024). ICLR 2024.
10. Love, M.I., et al. (2014). Genome Biology, 15, 550.
11. Mölder, F., et al. (2021). F1000Research, 10, 33.
12. Patro, R., et al. (2017). Nature Methods, 14, 417–419.
13. Stark, R., et al. (2019). Nature Reviews Genetics, 20, 631–656.
14. Subramanian, A., et al. (2005). PNAS, 102, 15545–15550.

---

## Supplementary Materials

- **Table S1:** Full DEG lists for each benchmark
- **Table S2:** Complete GSEA results (all collections, all benchmarks)
- **Table S3:** Decision logs from each benchmark run
- **Figure S1:** PCA plots for all benchmarks
- **Figure S2:** Volcano plots for all benchmarks
- **Figure S3:** GSEA barplots for all benchmarks
- **Figure S4:** Cross-method Venn diagrams
