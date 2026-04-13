# Supplementary Materials

## ARIA: Adaptive Reasoning for Integrated Analysis

---

### Table S1. Decision Log — Airway Benchmark

| Step | Decision Point | Trigger | Action | Confidence |
|------|---------------|---------|--------|------------|
| 1 | DP1 | QC metrics: min 15.2M reads | PASS — all samples above threshold | 0.95 |
| 2 | DP3 | 4 cell lines × 2 conditions | Use paired model `~ cell + dex` | 0.90 |
| 3 | DP2 | 951 DEGs at \|LFC\|>1 | Standard strategy (ORA + GSEA) | 0.95 |
| 4 | DP8 | Analysis complete | Generate HTML report + figures | 1.00 |

### Table S2. Decision Log — SEQC Benchmark

| Step | Decision Point | Trigger | Action | Confidence |
|------|---------------|---------|--------|------------|
| 1 | DP1 | QC metrics: min 37.8M counts | PASS | 0.95 |
| 2 | DP2 | 10,430 DEGs at \|LFC\|>1 | Standard strategy (EASY case) | 0.99 |
| 3 | DP8 | Analysis complete | Generate report | 1.00 |

### Table S3. Decision Log — Fmr1 KO Benchmark

| Step | Decision Point | Trigger | Action | Confidence |
|------|---------------|---------|--------|------------|
| 1 | DP1 | QC metrics: min 30.3M counts | PASS | 0.95 |
| 2 | DP2 | 398 DEGs at \|LFC\|>1 | Standard strategy | 0.90 |
| 3 | DP5 | Primary DE complete | Cross-validate: edgeR + limma-voom | 0.85 |
| 4 | DP8 | Analysis complete | Generate report | 1.00 |

### Table S4. Decision Log — Pasilla Benchmark

| Step | Decision Point | Trigger | Action | Confidence |
|------|---------------|---------|--------|------------|
| 1 | DP1 | QC metrics: min 8.4M counts | PASS | 0.95 |
| 2 | DP3 | Mixed SE/PE libraries detected | Include library type as covariate `~ type + condition` | 0.85 |
| 3 | DP2 | 224 DEGs at \|LFC\|>1 | Standard strategy | 0.90 |
| 4 | DP5 | Primary DE complete | Cross-validate: edgeR + limma-voom (r=0.9906) | 0.85 |
| 5 | DP7 | Results available | Sensitivity analysis: 12 cutoff combinations | 0.90 |
| 6 | DP8 | Analysis complete | Generate report | 1.00 |

### Table S5. Sensitivity Analysis — Pasilla

| padj cutoff | LFC cutoff | Number of DEGs |
|-------------|-----------|----------------|
| 0.01 | 2.0 | ~50 |
| 0.01 | 1.0 | ~150 |
| 0.01 | 0.5 | ~400 |
| 0.01 | 0.0 | ~700 |
| 0.05 | 2.0 | ~80 |
| 0.05 | 1.0 | 224 |
| 0.05 | 0.5 | 635 |
| 0.05 | 0.0 | 1,108 |
| 0.10 | 2.0 | ~100 |
| 0.10 | 1.0 | ~280 |
| 0.10 | 0.5 | ~800 |
| 0.10 | 0.0 | ~1,400 |

### Table S6. Known Target Gene Recovery

**Airway — Dexamethasone targets (7/7 = 100%)**

| Gene | Function | log2FC | padj | Detected |
|------|----------|--------|------|----------|
| DUSP1 | MAPK phosphatase | 2.94 | 1.3e-124 | Yes |
| KLF15 | Transcription factor | 4.46 | 7.8e-78 | Yes |
| CRISPLD2 | Anti-inflammatory | 2.63 | 1.6e-45 | Yes |
| PER1 | Circadian clock | 3.19 | 1.9e-81 | Yes |
| FKBP5 | GR co-chaperone | 4.04 | 7.4e-26 | Yes |
| TSC22D3 | Anti-inflammatory | 3.19 | 3.2e-19 | Yes |
| ZBTB16 | Transcription factor | 7.35 | 2.7e-40 | Yes |

**SEQC — Brain markers in HBRR (9/9 = 100%)**

| Gene | Function | log2FC (B vs A) | Correct direction | Detected |
|------|----------|----------------|-------------------|----------|
| GFAP | Astrocyte marker | +12.68 | Yes (UP in brain) | Yes |
| MBP | Oligodendrocyte | +8.00 | Yes | Yes |
| SYN1 | Synaptic | +7.54 | Yes | Yes |
| GAD1 | GABAergic | +4.35 | Yes | Yes |
| SLC17A7 | Glutamatergic | +8.14 | Yes | Yes |
| NEFL | Neurofilament | +7.20 | Yes | Yes |
| NEFM | Neurofilament | +8.14 | Yes | Yes |
| NEFH | Neurofilament | +3.78 | Yes | Yes |
| ENO2 | Neuron-specific | +2.19 | Yes | Yes |

**SEQC — Cancer markers in UHRR (6/7 = 86%)**

| Gene | log2FC (B vs A) | Correct direction | Detected |
|------|----------------|-------------------|----------|
| MYC | -5.69 | Yes (DOWN = UP in cancer) | Yes |
| CCND1 | -3.18 | Yes | Yes |
| CDK4 | -2.92 | Yes | Yes |
| EGFR | -0.44 | Yes | Yes |
| ERBB2 | -1.90 | Yes | Yes |
| TP53 | -3.39 | Yes | Yes |
| BCL2 | +0.37 | No (known brain expression) | Unexpected |

**Fmr1 KO — FMRP targets (9/11 detected)**

| Gene | Function | log2FC | padj | Detected |
|------|----------|--------|------|----------|
| Fmr1 | KO target gene | -0.983 | 3.7e-83 | Yes |
| Map1b | Microtubule | -0.573 | 7.1e-79 | Yes |
| Dlg4 | PSD-95 | -0.361 | 3.8e-22 | Yes |
| Shank3 | PSD scaffold | -0.454 | 2.3e-14 | Yes |
| Camk2a | Kinase | -0.424 | 7.4e-18 | Yes |
| Grin1 | NMDA-R | -0.919 | 2.5e-81 | Yes |
| Grin2a | NMDA-R | -0.812 | 1.1e-15 | Yes |
| Grin2b | NMDA-R | -0.479 | 1.0e-30 | Yes |
| Bdnf | Neurotrophin | -0.446 | 9.6e-03 | Yes |

### Table S7. Cross-Method Concordance

| Benchmark | DESeq2 DEGs | edgeR DEGs | limma DEGs | LFC r |
|-----------|------------|-----------|-----------|-------|
| SEQC | 13,818 | — | — | — |
| Airway | 2,426 | — | — | — |
| Fmr1 KO | 1,654 | ~1,600 | ~1,500 | 0.9999 |
| Pasilla | 635 | 739 | 674 | 0.9906 |

(DEG counts at padj<0.05, |LFC|>0.5)

### Figure S1. PCA plots for all benchmarks

Available in `benchmarks/results/[dataset]/DE/PCA.png`

### Figure S2. Volcano plots for all benchmarks

Available in `benchmarks/results/[dataset]/DE/Volcano.png`

### Figure S3. GSEA barplots for all benchmarks

Available in `benchmarks/results/[dataset]/GSEA/GSEA_*.png`

### Figure S4. Cross-method Venn diagrams

Available in `benchmarks/results/[dataset]/method_comparison/Venn_methods.png`

---

## Software Versions

| Component | Version |
|-----------|---------|
| Python | 3.10+ |
| R | 4.4.1 |
| DESeq2 | 1.46.0 |
| edgeR | 4.4.0 |
| limma | 3.62.1 |
| fgsea | 1.32.0 |
| msigdbr | 2026.1 |
| EnhancedVolcano | 1.24.0 |
| ggplot2 | 4.0.2 |
| Docker | 24.0+ |
| ARIA | 0.1.0 |
