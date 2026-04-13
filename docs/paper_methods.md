# 2. Methods

## 2.1 ARIA Architecture

ARIA is organized into four layers (Figure 1):

**Layer 1 — Reasoning Engine:** An LLM (Claude, Anthropic) serves as the reasoning backbone. It receives structured summaries of intermediate results (QC metrics, DEG counts, pathway enrichment statistics) and produces decisions in a standardized format specifying the action to take and the rationale.

**Layer 2 — Decision Protocol:** Eight formalized Decision Points (DPs) govern the analysis workflow. Each DP combines rule-based thresholds (e.g., mapping rate ≥ 85%) with LLM-driven contextual reasoning for edge cases. The protocol is implemented as a Python class (`ReasoningEngine`) with methods corresponding to each DP.

**Layer 3 — Execution Engine:** R and Python scripts are generated from parameterized templates and executed in Docker containers for reproducibility. The engine handles file I/O, Docker management, external API calls (STRING DB, MSigDB), and result parsing.

**Layer 4 — Analysis Modules:** Six validated bioinformatics modules:
- **DE Analysis:** DESeq2 (primary), with edgeR and limma-voom for cross-validation
- **GSEA:** fgsea with MSigDB collections (Hallmark, GO:BP/CC/MF, KEGG, Reactome)
- **ssGSEA:** Single-sample pathway activity scoring
- **PPI Network:** STRING DB protein-protein interaction analysis
- **Cell Type Deconvolution:** Marker-based cell composition estimation
- **WGCNA:** Weighted gene co-expression network analysis

## 2.2 Decision Protocol

### DP1: Quality Control Assessment
**Trigger:** Raw QC metrics available (from MultiQC or direct parsing)
**Rules:**
- Mapping rate < 85% → WARNING
- rRNA contamination > 5% → WARNING
- 5'→3' bias outside [0.8, 1.2] → WARNING
- Any sample flagged → LLM evaluates overall severity

**Action:** PASS (continue), WARNING (flag but continue), or FAIL (abort with explanation)

### DP2: DE Strategy Selection
**Trigger:** DESeq2 results available with DEG counts
**Rules:**
- DEG ≥ 100 (|LFC|>1) → Standard: ORA + GSEA
- DEG 10-100 → Add GSEA, try relaxed LFC (0.5), consider multi-factor model
- DEG < 10 → Prioritize GSEA; ORA likely underpowered

**Rationale:** Based on Subramanian et al. (2005), GSEA leverages the full ranked gene list and is more sensitive when individual gene effects are small but coordinated across pathways.

### DP3: Experimental Design Recognition
**Trigger:** Sample metadata available
**Rules:**
- Detect paired/blocked designs from sample naming or metadata
- Detect multi-factor designs (tissue + treatment)
- Select appropriate DESeq2 model formula

**Validation:** Airway benchmark demonstrates 21-47% DEG increase with correct paired model.

### DP4: Unexpected Signature Detection
**Trigger:** DE results contain gene clusters inconsistent with expected cell types
**Rules:**
- Ependymal markers (Foxj1, Dnah family) in neuronal tissue → flag, add deconvolution
- Immune markers in non-immune tissue → flag
- LLM evaluates whether pattern is biological or technical artifact

### DP5: Cross-Method Validation
**Trigger:** Primary DE analysis complete
**Action:** Run edgeR (exact test + QLF) and limma-voom on same data
**Output:** Concordance metrics (LFC correlation, DEG overlap Venn diagram)

### DP6: Literature-Based Interpretation
**Trigger:** Significant DEGs and pathways identified
**Action:** LLM searches its knowledge base for relevant prior studies
**Output:** Gene function annotations, pathway context, testable hypotheses
**Caveat:** Interpretations are explicitly labeled as hypothesis-level

### DP7: Sensitivity Analysis
**Trigger:** DEG results available
**Action:** Test multiple LFC cutoffs (1.0, 0.5, 0.0) and padj thresholds (0.05, 0.1)
**Output:** Sensitivity table showing DEG counts across thresholds

### DP8: Report Generation
**Trigger:** All analyses complete
**Action:** Generate HTML report with embedded figures, Excel data files, PNG figures
**Output:** Publication-ready report adapted to specified audience level

## 2.3 Implementation Details

ARIA is implemented in Python 3.10+ with analysis modules executed as R scripts within Docker containers. Key dependencies:

| Component | Tool | Version |
|-----------|------|---------|
| DE Analysis | DESeq2 | 1.46.0 |
| DE Validation | edgeR, limma | 4.0+, 3.62+ |
| GSEA | fgsea | 1.32.0 |
| Gene Sets | msigdbr | 2026.1 |
| Network | STRING DB API | v12 |
| Co-expression | WGCNA | 1.73 |
| LFC Shrinkage | apeglm | 1.28.0 |
| P-value Adjustment | IHW | 1.34.0 |
| Annotation | org.Mm.eg.db, org.Hs.eg.db | 3.20.0 |

All R scripts are generated from parameterized templates (`aria/modules/`), ensuring reproducibility and enabling parameter sweeps.

## 2.4 Benchmark Datasets

| Dataset | GEO | Species | Design | N | Expected DEGs | Tests |
|---------|-----|---------|--------|---|---------------|-------|
| SEQC/MAQC-III | GSE49712 | Human | A vs B reference RNA | 10 | ~8,000 | DP1, DP2 (easy) |
| Bottomly | GSE26024 | Mouse | C57BL/6J vs DBA/2J brain | 21 | ~300 | DP2 (moderate), DP5 |
| Airway | GSE52778 | Human | Dex ± ASM cells, paired | 8 | ~2,000 | DP1, DP2, DP3 |
| Fmr1 KO | GSE84989 | Mouse | KO vs WT brain | ~8 | ~100 | DP2 (hard), DP4, DP6 |

## 2.5 Evaluation Metrics

**A. Analysis Accuracy:**
- DEG Precision and Recall vs. ground truth (SEQC qRT-PCR) or expert consensus
- Known target gene recovery rate (Airway: 7 dexamethasone targets)
- log2FC Pearson correlation across methods

**B. Decision Quality:**
- Strategy Appropriateness: Whether DP decisions match expert recommendations
- Gain from Correct Decisions: Quantified as % increase in DEGs (e.g., paired vs unpaired)

**C. Practical Utility:**
- Analysis wall-clock time
- Number of output files and report completeness
- API cost per analysis

**D. Reproducibility:**
- Triplicate execution concordance
- Docker-based full reproducibility
