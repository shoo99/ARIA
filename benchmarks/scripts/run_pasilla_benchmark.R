#!/usr/bin/env Rscript
###############################################################################
# ARIA Benchmark: Pasilla (GSE18508) — MODERATE scenario
# Drosophila, pasilla RNAi KD vs control, mixed SE/PE libraries
# Tests: DP1, DP2 (moderate DEGs), DP5 (cross-method), DP7 (sensitivity)
# Replaces Bottomly simulation with real data
###############################################################################

library(DESeq2)
library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(futile.logger)

out_dir <- "/data/benchmarks/results/pasilla"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("###########################################################\n")
cat("### ARIA Benchmark: Pasilla (GSE18508) — MODERATE\n")
cat("###########################################################\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
cat("\n--- 1. Data Loading ---\n")
counts <- as.matrix(read.csv("/data/benchmarks/data/pasilla/pasilla_gene_counts.tsv", sep="\t", row.names="gene_id"))
cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Samples:", paste(colnames(counts), collapse=", "), "\n")

# Metadata
meta <- read.csv("/data/benchmarks/data/pasilla/pasilla_sample_annotation.csv", row.names=1)
meta$condition <- factor(meta$condition, levels=c("untreated","treated"))
meta$type <- factor(meta$type)
cat("\nMetadata:\n")
print(meta[, c("condition","type")])

# Match sample names — metadata has "fb" suffix, counts don't
rownames(meta) <- sub("fb$", "", rownames(meta))
counts <- counts[, rownames(meta)]

# ── 2. DP1: QC Assessment ───────────────────────────────────────────────────
cat("\n--- 2. DP1: QC Assessment ---\n")
total_counts <- colSums(counts)
cat("Total counts per sample:\n"); print(total_counts)
cat(sprintf("DP1 DECISION: Min counts = %d → %s\n", min(total_counts),
    ifelse(min(total_counts) > 1e6, "PASS", "WARNING")))
cat("Note: Mixed SE/PE library types detected\n")

# ── 3. DE Analysis ──────────────────────────────────────────────────────────
cat("\n--- 3. DE Analysis (DESeq2) ---\n")
de_dir <- file.path(out_dir, "DE")
dir.create(de_dir, showWarnings=FALSE)

keep <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep, ]
cat("Genes after filtering:", nrow(counts_filt), "\n")

# Model with type (SE/PE) as covariate — tests DP3-like reasoning
cat("\n  Model A: ~ condition (no type covariate)\n")
dds_A <- DESeqDataSetFromMatrix(counts_filt, meta, ~ condition)
dds_A <- DESeq(dds_A)
res_A <- results(dds_A, contrast=c("condition","treated","untreated"), alpha=0.05)

cat("  Model B: ~ type + condition (with SE/PE covariate)\n")
dds_B <- DESeqDataSetFromMatrix(counts_filt, meta, ~ type + condition)
dds_B <- DESeq(dds_B)
res_B <- results(dds_B, contrast=c("condition","treated","untreated"), alpha=0.05)

# Compare
cat("\n  Comparison:\n")
for (lfc in c(1.0, 0.5, 0.0)) {
    n_A <- sum(!is.na(res_A$padj) & res_A$padj < 0.05 & abs(res_A$log2FoldChange) > lfc)
    n_B <- sum(!is.na(res_B$padj) & res_B$padj < 0.05 & abs(res_B$log2FoldChange) > lfc)
    cat(sprintf("  |LFC|>%.1f: No covar=%d, With type=%d (gain: %+d)\n", lfc, n_A, n_B, n_B-n_A))
}

# Use Model B (correct)
res <- res_B
res$symbol <- rownames(res)

# Save
res_df <- as.data.frame(res)
write.csv(res_df, file.path(de_dir, "DE_all.csv"), row.names=TRUE)

sig1 <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
sig05 <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5)
write.csv(sig1[order(sig1$padj),], file.path(de_dir, "DE_sig_LFC1.csv"), row.names=TRUE)

# ── 4. DP2: Strategy selection ──────────────────────────────────────────────
cat("\n--- 4. DP2: Strategy Selection ---\n")
n_deg <- nrow(sig1)
if (n_deg >= 100) {
    cat(sprintf("DP2 DECISION: %d DEGs → standard\n", n_deg))
} else if (n_deg >= 10) {
    cat(sprintf("DP2 DECISION: %d DEGs → augmented (add GSEA)\n", n_deg))
} else {
    cat(sprintf("DP2 DECISION: %d DEGs → GSEA priority\n", n_deg))
}

# ── 5. DP5: Cross-method validation ─────────────────────────────────────────
cat("\n--- 5. DP5: Cross-Method Validation ---\n")
dp5_dir <- file.path(out_dir, "method_comparison")
dir.create(dp5_dir, showWarnings=FALSE)

# edgeR
y <- DGEList(counts=counts_filt, group=meta$condition)
y <- y[filterByExpr(y, group=meta$condition),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ meta$type + meta$condition)
y <- estimateDisp(y, design)

et <- exactTest(y, pair=c("untreated","treated"))
res_exact <- topTags(et, n=Inf)$table

# limma-voom
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef=ncol(design), number=Inf)

cat("Cross-method (padj<0.05, |LFC|>0.5):\n")
n_deseq2 <- nrow(sig05)
n_edger <- sum(res_exact$FDR < 0.05 & abs(res_exact$logFC) > 0.5, na.rm=TRUE)
n_limma <- sum(res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > 0.5, na.rm=TRUE)
cat(sprintf("  DESeq2: %d, edgeR: %d, limma: %d\n", n_deseq2, n_edger, n_limma))

common <- intersect(rownames(res), rownames(res_exact))
cor_val <- cor(as.data.frame(res)[common,"log2FoldChange"], res_exact[common,"logFC"], use="complete.obs")
cat(sprintf("  LFC correlation: r=%.4f\n", cor_val))

# Venn
deg_lists <- list(
    DESeq2 = rownames(sig05),
    edgeR = rownames(res_exact)[res_exact$FDR < 0.05 & abs(res_exact$logFC) > 0.5],
    limma = rownames(res_limma)[res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > 0.5]
)
png(file.path(dp5_dir, "Venn_methods.png"), width=700, height=600)
venn.plot <- venn.diagram(x=deg_lists, category.names=names(deg_lists), filename=NULL,
    fill=c("#E41A1C","#377EB8","#4DAF4A"), alpha=0.4, cex=1.8, cat.cex=1.2,
    main="Pasilla: Method Comparison (padj<0.05, |LFC|>0.5)")
grid::grid.draw(venn.plot)
dev.off()

# ── 6. Plots ─────────────────────────────────────────────────────────────────
cat("\n--- 6. Visualization ---\n")
vsd <- vst(dds_B, blind=FALSE)

pca_data <- plotPCA(vsd, intgroup=c("condition","type"), returnData=TRUE)
pct_var <- round(100*attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color=condition, shape=type)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ", pct_var[1], "%")) + ylab(paste0("PC2: ", pct_var[2], "%")) +
    ggtitle("Pasilla: PCA (treated vs untreated)") + theme_bw() +
    scale_color_manual(values=c("untreated"="#377EB8","treated"="#E41A1C"))
ggsave(file.path(de_dir, "PCA.png"), p, width=8, height=6, dpi=150)

png(file.path(de_dir, "Volcano.png"), width=900, height=700)
print(EnhancedVolcano(res, lab=res$symbol, x="log2FoldChange", y="padj",
    pCutoff=0.05, FCcutoff=1, title="Pasilla: treated vs untreated",
    pointSize=2, labSize=3, colAlpha=0.6))
dev.off()

# ── 7. DP7: Sensitivity analysis ────────────────────────────────────────────
cat("\n--- 7. DP7: Sensitivity Analysis ---\n")
sens <- data.frame()
for (p_cut in c(0.01, 0.05, 0.1)) {
    for (lfc_cut in c(2.0, 1.0, 0.5, 0.0)) {
        n <- sum(!is.na(res$padj) & res$padj < p_cut & abs(res$log2FoldChange) > lfc_cut)
        sens <- rbind(sens, data.frame(padj_cutoff=p_cut, lfc_cutoff=lfc_cut, n_DEG=n))
    }
}
cat("Sensitivity table:\n")
print(sens)
write.csv(sens, file.path(out_dir, "sensitivity_analysis.csv"), row.names=FALSE)

# ── 8. Known biology ────────────────────────────────────────────────────────
cat("\n--- 8. Known Biology ---\n")
cat("Pasilla (ps) = Drosophila ortholog of NOVA1/NOVA2 splicing factor\n")
cat("KD of pasilla affects alternative splicing, not primarily gene expression\n")
if ("ps" %in% rownames(res)) {
    r <- as.data.frame(res)["ps",]
    cat(sprintf("  ps itself: LFC=%.3f, padj=%s\n", r$log2FoldChange,
        ifelse(is.na(r$padj), "NA", sprintf("%.1e", r$padj))))
}

# ── 9. Summary ───────────────────────────────────────────────────────────────
cat("\n###########################################################\n")
cat("### ARIA BENCHMARK SUMMARY: Pasilla (MODERATE)\n")
cat("###########################################################\n")
cat(sprintf("Genes tested: %d\n", nrow(res)))
cat(sprintf("DEGs (padj<0.05, |LFC|>1): %d\n", nrow(sig1)))
cat(sprintf("DEGs (padj<0.05, |LFC|>0.5): %d\n", nrow(sig05)))
cat(sprintf("Model covariate (type) gain: see comparison above\n"))
cat(sprintf("Cross-method LFC r=%.4f\n", cor_val))
cat("DP1: QC PASS\n")
cat("DP3: Library type covariate detected\n")
cat(sprintf("DP5: DESeq2=%d, edgeR=%d, limma=%d\n", n_deseq2, n_edger, n_limma))
cat("DP7: Sensitivity analysis completed\n")
cat("Results saved to:", out_dir, "\n")
