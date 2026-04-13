# bioRxiv Submission Checklist

## Manuscript Information

- **Title:** ARIA: Adaptive Reasoning for Integrated Analysis — An LLM-Powered Framework for Autonomous Transcriptome Analysis with Decision-Aware Workflow Orchestration
- **Subject Area:** Bioinformatics
- **Article Type:** New Results

## Files to Submit

### Required
1. **Manuscript PDF** — Convert `manuscript_full.md` to PDF (via Pandoc or LaTeX)
2. **Figures:**
   - Figure 1: Architecture (`docs/figures/Fig1_architecture.pdf`)
   - Figure 2: Decision Flow (`docs/figures/Fig2_decision_flow.pdf`)
   - Figure 3: Airway Benchmark (`docs/figures/Fig3_airway_benchmark.pdf`)
   - Figure 4: Comparison Table (`docs/figures/Fig4_comparison.pdf`)

### Supplementary
3. **Supplementary Materials** — `supplementary.md` → PDF

## Conversion Commands

```bash
# Install pandoc if needed
# apt install pandoc texlive-latex-recommended texlive-fonts-recommended

# Convert manuscript to PDF
pandoc docs/manuscript_full.md -o docs/ARIA_manuscript.pdf \
  --pdf-engine=xelatex \
  -V geometry:margin=1in \
  -V fontsize=11pt \
  -V mainfont="Times New Roman" \
  --number-sections

# Convert supplementary to PDF
pandoc docs/supplementary.md -o docs/ARIA_supplementary.pdf \
  --pdf-engine=xelatex \
  -V geometry:margin=1in \
  -V fontsize=10pt \
  --number-sections
```

## bioRxiv Posting Criteria

✅ Full research paper with methodological details — YES (4 layers, 8 DPs, 6 modules)
✅ Contains results — YES (4 benchmarks with quantitative metrics)
✅ Not a simple automated analysis of public data — YES (novel framework + systematic evaluation)
✅ Includes methods — YES (architecture, implementation, evaluation metrics)
✅ Software availability — YES (GitHub + Docker)

## Pre-submission Checks

- [ ] All figures are high resolution (300 dpi minimum)
- [ ] All benchmark data sources are cited with GEO accessions
- [ ] GitHub repository is public and accessible
- [ ] Docker image is available
- [ ] License is specified (MIT)
- [ ] No proprietary/client data is included
- [ ] All references are complete
- [ ] Abstract is under 300 words
- [ ] Author information is complete
- [ ] Conflict of interest statement is included

## Notes

- bioRxiv does not require peer review
- Once posted, the preprint gets a DOI immediately
- Can be updated with new versions
- Can submit to a journal simultaneously
- Recommended journals: Bioinformatics (Oxford), NAR (Software), GigaScience
