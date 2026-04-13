# ARIA: Adaptive Reasoning for Integrated Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**An LLM-powered autonomous framework for transcriptome analysis with decision-aware workflow orchestration.**

ARIA goes beyond traditional fixed pipelines by using a Large Language Model as a reasoning engine that evaluates intermediate results and adaptively selects the next analysis step — mimicking how an experienced bioinformatician thinks.

## Key Features

- **Adaptive Decision Making**: Automatically adjusts analysis strategy based on intermediate results (e.g., switches to GSEA when DEGs are insufficient)
- **Multi-method Validation**: Cross-validates results using DESeq2, edgeR, and limma-voom
- **Comprehensive Reporting**: Generates publication-ready HTML reports with embedded visualizations
- **Reproducible**: Docker-based execution ensures full reproducibility
- **Extensible**: Modular architecture allows adding new analysis tools

## Architecture

```
┌─────────────────────────────────────────┐
│  Layer 1: Reasoning Engine (LLM)        │
│  - QC evaluation & decision making      │
│  - Strategy adaptation                  │
│  - Literature-based interpretation      │
├─────────────────────────────────────────┤
│  Layer 2: Decision Protocol             │
│  - 8 Decision Points (DP1-DP8)          │
│  - Rule-based + LLM-reasoning hybrid    │
├─────────────────────────────────────────┤
│  Layer 3: Execution Engine              │
│  - Docker container management          │
│  - R/Python script generation & run     │
│  - External API calls                   │
├─────────────────────────────────────────┤
│  Layer 4: Analysis Modules              │
│  - DESeq2, edgeR, limma-voom           │
│  - fgsea, WGCNA, STRING DB             │
│  - Cell type deconvolution              │
└─────────────────────────────────────────┘
```

## Decision Points

| DP | Trigger | Action |
|----|---------|--------|
| DP1 | QC metrics available | Pass/Fail/Warning assessment |
| DP2 | DEG count < threshold | Switch to GSEA, adjust LFC cutoff |
| DP3 | Multi-condition design | Choose per-group or multi-factor model |
| DP4 | Unexpected signature | Add cell type deconvolution |
| DP5 | Method disagreement | Cross-validate with multiple DE tools |
| DP6 | Pathway results | Literature search + hypothesis generation |
| DP7 | Parameter uncertainty | Sensitivity analysis |
| DP8 | Analysis complete | Generate audience-appropriate report |

## Quick Start

```bash
# Install
pip install aria-omics

# Run with default settings
aria analyze --input samplesheet.csv --genome GRCm39 --output results/

# Run with specific LLM backend
aria analyze --input samplesheet.csv --llm claude --output results/
```

## Docker

```bash
docker pull shoo99/aria:latest
docker run -v $(pwd):/data aria:latest analyze --input /data/samplesheet.csv
```

## Benchmarks

ARIA has been evaluated on 4 public RNA-seq datasets spanning different difficulty levels:

| Dataset | Accession | Difficulty | DEGs | Use Case |
|---------|-----------|------------|------|----------|
| SEQC | GSE49712 | Easy | ~8,000 | Accuracy (qRT-PCR ground truth) |
| Bottomly | GSE26024 | Moderate | ~300 | Reproducibility (mouse brain) |
| Airway | GSE52778 | Complex | ~3,000 | Design recognition (paired) |
| Fmr1 KO | GSE84989 | Hard | ~100 | Adaptive strategy switching |

## Citation

If you use ARIA in your research, please cite:

```bibtex
@article{kang2026aria,
  title={ARIA: Adaptive Reasoning for Integrated Analysis — An LLM-Powered Framework for Autonomous Transcriptome Analysis with Decision-Aware Workflow Orchestration},
  author={Kang, Byeongsoo},
  journal={bioRxiv},
  year={2026},
  doi={pending}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
