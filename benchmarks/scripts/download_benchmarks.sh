#!/bin/bash
# Download benchmark datasets for ARIA evaluation
# All datasets are publicly available on GEO/SRA

set -euo pipefail

DATA_DIR="$(dirname "$0")/../data"
mkdir -p "$DATA_DIR"

echo "=== ARIA Benchmark Data Download ==="
echo "This script downloads public RNA-seq datasets for benchmarking."
echo ""

# 1. SEQC/MAQC-III (GSE49712) - EASY
echo "--- [1/4] SEQC (GSE49712) ---"
echo "Human reference RNA, 5 Sample A vs 5 Sample B"
echo "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49712"
# TODO: Add SRA download commands (prefetch + fasterq-dump)

# 2. Bottomly (GSE26024) - MODERATE
echo "--- [2/4] Bottomly (GSE26024) ---"
echo "Mouse brain, C57BL/6J vs DBA/2J, 10 vs 11 samples"
echo "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26024"

# 3. Airway (GSE52778) - MODERATE + COMPLEX
echo "--- [3/4] Airway (GSE52778) ---"
echo "Human ASM cells, Dex-treated vs control, paired design"
echo "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778"
echo "Also available as Bioconductor 'airway' package"

# 4. Fmr1 KO (GSE84989) - HARD
echo "--- [4/4] Fmr1 KO (GSE84989) ---"
echo "Mouse brain, Fmr1 KO vs WT"
echo "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84989"

echo ""
echo "=== Download instructions ==="
echo "Use SRA Toolkit to download FASTQ files:"
echo "  prefetch SRRXXXXXXX"
echo "  fasterq-dump SRRXXXXXXX --outdir \$DATA_DIR"
echo ""
echo "Or use nf-core/fetchngs pipeline:"
echo "  nextflow run nf-core/fetchngs --input ids.csv --outdir \$DATA_DIR"
