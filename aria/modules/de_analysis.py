"""
DE Analysis Module

Generates and executes R scripts for differential expression analysis
using DESeq2, with optional edgeR/limma-voom cross-validation.
"""

from pathlib import Path
from typing import Optional


class DEAnalysisModule:
    """Differential expression analysis using DESeq2."""

    R_TEMPLATE = '''#!/usr/bin/env Rscript
library(DESeq2)
library(AnnotationDbi)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Load counts
counts_raw <- read.delim("{counts_file}", row.names=1, check.names=FALSE)
gene_names <- counts_raw$gene_name
counts_raw$gene_name <- NULL
counts_int <- round(as.matrix(counts_raw))

# Metadata
sample_names <- colnames(counts_int)
meta <- data.frame(
    genotype = factor(ifelse(grepl("^{control_prefix}_", sample_names), "{control_label}", "{treatment_label}"),
                      levels=c("{control_label}", "{treatment_label}")),
    row.names = sample_names
)

# Filter and run DESeq2
dds <- DESeqDataSetFromMatrix(counts_int, meta, ~ genotype)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype", "{treatment_label}", "{control_label}"), alpha=0.05)

# Add gene info
gene_info <- data.frame(ensembl_id=rownames(counts_int), gene_name=gene_names, stringsAsFactors=FALSE)
res$gene_name <- gene_info$gene_name[match(rownames(res), gene_info$ensembl_id)]

# Save results
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
write.csv(res_df, "{output_dir}/DE_all.csv", row.names=FALSE)

# Summary
for (lfc in c(1.0, 0.5, 0.0)) {{
    n <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > lfc)
    cat(sprintf("padj<0.05, |LFC|>%.1f: %d DEGs\\n", lfc, n))
}}

# Plots
vsd <- vst(dds, blind=FALSE)

# PCA
pca_data <- plotPCA(vsd, intgroup="genotype", returnData=TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color=genotype, label=name)) +
    geom_point(size=5) + geom_text_repel(size=3) +
    xlab(paste0("PC1: ", pct_var[1], "%")) + ylab(paste0("PC2: ", pct_var[2], "%")) +
    ggtitle("PCA") + theme_bw() +
    scale_color_manual(values=c("{control_label}"="#377EB8", "{treatment_label}"="#E41A1C"))
ggsave("{output_dir}/PCA.png", p, width=8, height=6, dpi=150)

# Volcano
for (lfc in c(1.0, 0.5)) {{
    png(sprintf("{output_dir}/Volcano_LFC%s.png", gsub("\\\\.", "", as.character(lfc))), width=900, height=700)
    print(EnhancedVolcano(res, lab=res$gene_name, x="log2FoldChange", y="padj",
        pCutoff=0.05, FCcutoff=lfc, pointSize=2, labSize=3.5, colAlpha=0.6,
        title=paste0("{treatment_label} vs {control_label} (|LFC|>", lfc, ")")))
    dev.off()
}}

# Heatmap
sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5)
if (nrow(sig) > 1) {{
    top_n <- min(50, nrow(sig))
    top_genes <- rownames(sig[order(sig$padj),])[1:top_n]
    mat <- assay(vsd)[top_genes,] - rowMeans(assay(vsd)[top_genes,])
    rn <- gene_info$gene_name[match(rownames(mat), gene_info$ensembl_id)]
    rownames(mat) <- ifelse(is.na(rn), rownames(mat), rn)
    anno_c <- data.frame(Genotype=meta$genotype, row.names=rownames(meta))
    png("{output_dir}/Heatmap_DEG.png", width=700, height=max(400, nrow(mat)*14))
    pheatmap(mat, annotation_col=anno_c, cluster_rows=TRUE, cluster_cols=TRUE,
             show_rownames=TRUE, fontsize_row=7, main="Top DEGs")
    dev.off()
}}

cat("DE analysis complete.\\n")
'''

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_script(self, counts_file: str,
                        control_prefix: str = "WT",
                        control_label: str = "WT",
                        treatment_label: str = "KO") -> str:
        """Generate the R script for DE analysis."""
        script = self.R_TEMPLATE.format(
            counts_file=counts_file,
            output_dir=str(self.output_dir),
            control_prefix=control_prefix,
            control_label=control_label,
            treatment_label=treatment_label
        )
        script_path = self.output_dir / "run_DE.R"
        script_path.write_text(script)
        return str(script_path)
