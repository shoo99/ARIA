# Research Square Submission Checklist — ARIA

**Submission URL:** https://www.researchsquare.com/browse/user-submit
**Estimated time:** 30-45 minutes
**Expected posting:** 1-2 business days after internal screening

---

## 1. Pre-filled Metadata (copy-paste ready)

### Title
```
ARIA: Adaptive Reasoning for Integrated Analysis — An LLM-Powered Framework for Autonomous Transcriptome Analysis with Decision-Aware Workflow Orchestration
```

### Short Title (for header, ≤70 chars)
```
ARIA: LLM-Powered Autonomous RNA-seq Analysis
```

### Article Type
`Research Article` (or `Methods Article` if offered)

### Subject Categories (select primary + secondaries)
- **Primary:** Bioinformatics
- **Secondary 1:** Computational Biology
- **Secondary 2:** Artificial Intelligence
- **Secondary 3:** Genomics

### Keywords (6)
```
RNA-seq, differential expression, large language models, autonomous analysis, decision-aware pipeline, transcriptomics
```

### Corresponding Author
- **Name:** Byeongsoo Kang (강병수)
- **Email:** shoo99@gmail.com
- **Affiliation:** SYSOFT, Republic of Korea
- **ORCID:** `0009-0007-2324-2351` (https://orcid.org/0009-0007-2324-2351)

### License
`CC BY 4.0` (권장 — 가장 유연, Research Square 기본값)

---

## 2. Abstract (paste into submission form)

*Word count: ~290 / 300 limit*

RNA-seq transcriptome analysis requires a multi-step workflow involving quality control, alignment, quantification, differential expression testing, pathway analysis, and biological interpretation. While automated pipelines such as nf-core/rnaseq execute these steps reproducibly, the critical decisions between steps — evaluating quality metrics, selecting statistical methods, adapting analysis strategies based on intermediate results, and interpreting findings in biological context — remain dependent on expert bioinformaticians. Here we present ARIA (Adaptive Reasoning for Integrated Analysis), an open-source framework that uses a Large Language Model (LLM) as a reasoning engine to autonomously navigate the decision space of RNA-seq analysis. ARIA implements eight Decision Points (DPs) that govern quality assessment, strategy adaptation, method selection, and result interpretation, combining rule-based thresholds with LLM-driven contextual reasoning. We benchmark ARIA on four public RNA-seq datasets spanning three species: SEQC (GSE49712, human, 10,430 DEGs), Airway (GSE52778, human paired design, 951 DEGs), Fmr1 KO (GSE180135, mouse, 398 DEGs), and Pasilla (GSE18508, Drosophila, 224 DEGs with mixed library types). ARIA correctly identifies paired experimental designs (increasing DEG detection by 21–47%), detects technical covariates such as library type (+4–30% DEG gain), adaptively selects analysis strategies based on DEG counts, and cross-validates results across DESeq2/edgeR/limma-voom (r > 0.99). All 7/7 known dexamethasone-responsive genes were recovered in Airway, 9/11 FMRP translational targets in Fmr1 KO, and 15/16 tissue-type markers in SEQC were correctly assigned. ARIA is freely available at https://github.com/shoo99/ARIA under the MIT license.

---

## 3. Files to Upload

| Required | File | Size | Path |
|---|---|---|---|
| ✅ Main manuscript PDF | `ARIA_bioRxiv.pdf` | 1.2 MB | `docs/ARIA_bioRxiv.pdf` |
| ✅ Supplementary PDF | `ARIA_supplementary.pdf` | 28 KB | `docs/ARIA_supplementary.pdf` |
| Optional | Korean version | 1.8 MB | `docs/ARIA_korean.pdf` |
| Optional | Figures (high-res) | — | `docs/figures/Fig1-4.pdf` |

*Research Square는 main PDF에 figure가 embedded되어 있으면 추가 figure 파일 업로드 불필요.*

---

## 4. Cover Letter (optional but recommended)

### Draft Cover Letter

```
Dear Research Square Editorial Team,

I am submitting the manuscript titled "ARIA: Adaptive Reasoning for Integrated
Analysis — An LLM-Powered Framework for Autonomous Transcriptome Analysis with
Decision-Aware Workflow Orchestration" for consideration as a preprint.

ARIA introduces a novel framework that uses a Large Language Model as a
reasoning engine to autonomously navigate the multi-step decision space of
RNA-seq transcriptome analysis. While current workflow managers (Nextflow,
Snakemake) and curated pipelines (nf-core/rnaseq) automate execution, the
critical decisions between steps — quality assessment, design recognition,
strategy adaptation, method selection, and biological interpretation — still
require expert bioinformaticians. ARIA addresses this gap by formalizing
eight Decision Points and combining rule-based thresholds with LLM-driven
contextual reasoning.

To our knowledge, this is the first systematic application of LLM-based
adaptive reasoning to the complete RNA-seq analysis workflow. We validate
ARIA on four public datasets spanning three species and diverse complexity
levels, demonstrating robust recovery of known biological markers and
adaptive behavior in response to data characteristics (paired designs,
library-type heterogeneity, low-DEG scenarios).

The framework is open-source (MIT license) with full source code available
at https://github.com/shoo99/ARIA, and is reproducible via Docker. Benchmark
datasets are all publicly available from GEO (accessions listed in Methods).

This preprint provides a foundation for community feedback prior to
peer-reviewed journal submission. There are no conflicts of interest to
declare. All data sources are public and properly cited.

Thank you for considering this submission.

Sincerely,
Byeongsoo Kang
SYSOFT, Republic of Korea
shoo99@gmail.com
```

---

## 5. Conflict of Interest Statement
```
The author declares no competing interests.
```

---

## 6. Funding Statement
```
No external funding was received for this work.
```
*(SYSOFT 내부 자금이면 해당 내용 추가)*

---

## 7. Data Availability Statement
```
All benchmark datasets are publicly available from the Gene Expression Omnibus
(GEO) under accession numbers GSE49712 (SEQC), GSE52778 (Airway), GSE180135
(Fmr1 KO), and GSE18508 (Pasilla). ARIA source code is available at
https://github.com/shoo99/ARIA under the MIT license. Docker images are
available for full computational reproducibility.
```

---

## 8. Author Contributions
```
B.K. conceived the framework, implemented the software, performed all
analyses, and wrote the manuscript.
```

---

## 9. Ethics Statement
```
This study exclusively used publicly available de-identified datasets from
the Gene Expression Omnibus. No new animal or human subjects research was
performed. Institutional ethics approval was therefore not required.
```

---

## 10. Pre-submission Checklist

- [ ] GitHub repo (shoo99/ARIA) public and accessible
- [ ] MIT LICENSE file present
- [ ] All benchmark GEO accessions verified and cited
- [ ] Abstract ≤300 words (confirmed: 290 words)
- [ ] Corresponding author email checked
- [x] ORCID acquired: 0009-0007-2324-2351
- [ ] Figures in PDF readable at 100% zoom
- [ ] No client/proprietary data referenced (NGS8 not mentioned)
- [ ] All reference citations complete
- [ ] PDF font embedded (check with `pdffonts ARIA_bioRxiv.pdf`)

---

## 11. Post-submission Next Steps

1. **1-2 days:** Editorial screening (auto-check of formatting, scope)
2. **Posting:** DOI assigned (`10.21203/rs.3.rs-XXXXXX/v1`)
3. **Citation format:** Kang, B. "ARIA: ..." Research Square (2026). doi:XXX
4. **Update README.md** in GitHub to include preprint DOI badge
5. **Parallel:** Continue JOSS submission prep (software paper)
6. **Later:** Submit to peer-reviewed journal (Bioinformatics Oxford, SoftwareX, or NAR)

---

## 12. Update GitHub README after posting

Add badge to top of README.md:
```markdown
[![Preprint](https://img.shields.io/badge/Research%20Square-10.21203%2Frs.3.rs--XXXXX-blue)](https://doi.org/10.21203/rs.3.rs-XXXXX/v1)
```

And citation block:
```markdown
## Citation

If you use ARIA in your research, please cite:

> Kang, B. (2026). ARIA: Adaptive Reasoning for Integrated Analysis — An
> LLM-Powered Framework for Autonomous Transcriptome Analysis with
> Decision-Aware Workflow Orchestration. *Research Square*.
> https://doi.org/10.21203/rs.3.rs-XXXXX/v1
```
