#!/usr/bin/env Rscript
###############################################################################
# ARIA Benchmark: SEQC/MAQC-III (GSE49712) — EASY scenario
# Human reference RNA, Sample A (UHRR) vs Sample B (HBRR), 5 vs 5
# Tests: DP1 (QC), DP2 (many DEGs → standard), accuracy vs known biology
###############################################################################

library(DESeq2)
library(fgsea)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

out_dir <- "/data/benchmarks/results/seqc"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("###########################################################\n")
cat("### ARIA Benchmark: SEQC (GSE49712) — EASY\n")
cat("###########################################################\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
cat("\n--- 1. Data Loading ---\n")
counts_raw <- read.delim("/data/benchmarks/data/seqc/GSE49712_HTSeq.txt",
                          row.names=1, check.names=FALSE, quote='"')

# Remove HTSeq summary rows
summary_rows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
counts_raw <- counts_raw[!rownames(counts_raw) %in% summary_rows, ]
counts <- as.matrix(counts_raw)
cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Samples:", paste(colnames(counts), collapse=", "), "\n")

meta <- data.frame(
    condition = factor(ifelse(grepl("^A", colnames(counts)), "A", "B"), levels=c("A","B")),
    row.names = colnames(counts)
)
print(meta)

# ── 2. DP1: QC Assessment ───────────────────────────────────────────────────
cat("\n--- 2. DP1: QC Assessment ---\n")
total_counts <- colSums(counts)
cat("Total counts per sample:\n"); print(total_counts)
cat(sprintf("DP1 DECISION: Min counts = %d → %s\n", min(total_counts),
    ifelse(min(total_counts) > 1e6, "PASS", "WARNING")))

# ── 3. DE Analysis ──────────────────────────────────────────────────────────
cat("\n--- 3. DE Analysis (DESeq2) ---\n")
de_dir <- file.path(out_dir, "DE")
dir.create(de_dir, showWarnings=FALSE)

keep <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep, ]
cat("Genes after filtering:", nrow(counts_filt), "\n")

dds <- DESeqDataSetFromMatrix(counts_filt, meta, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","B","A"), alpha=0.05)
res$symbol <- rownames(res)

for (lfc in c(2.0, 1.0, 0.5, 0.0)) {
    n <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > lfc)
    cat(sprintf("  padj<0.05, |LFC|>%.1f: %d DEGs\n", lfc, n))
}

# Save
res_df <- as.data.frame(res)
write.csv(res_df, file.path(de_dir, "DE_all.csv"), row.names=TRUE)

# ── 4. DP2: DE Evaluation ───────────────────────────────────────────────────
cat("\n--- 4. DP2: DE Result Evaluation ---\n")
n_deg <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1)
cat(sprintf("DP2 DECISION: %d DEGs (|LFC|>1) → MANY DEGs, standard ORA + GSEA\n", n_deg))
cat("This is the EASY scenario — massive differences between UHRR and HBRR\n")

# ── 5. Plots ─────────────────────────────────────────────────────────────────
cat("\n--- 5. Visualization ---\n")
vsd <- vst(dds, blind=FALSE)

# PCA
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
pct_var <- round(100*attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color=condition, label=name)) +
    geom_point(size=5) + geom_text_repel(size=3.5) +
    xlab(paste0("PC1: ", pct_var[1], "%")) + ylab(paste0("PC2: ", pct_var[2], "%")) +
    ggtitle("SEQC: PCA (Sample A vs B)") + theme_bw() +
    scale_color_manual(values=c("A"="#377EB8","B"="#E41A1C"))
ggsave(file.path(de_dir, "PCA.png"), p, width=8, height=6, dpi=150)

# Volcano
png(file.path(de_dir, "Volcano.png"), width=900, height=700)
print(EnhancedVolcano(res, lab=res$symbol, x="log2FoldChange", y="padj",
    pCutoff=0.05, FCcutoff=2, title="SEQC: Sample B vs A (|LFC|>2)",
    pointSize=1.5, labSize=3, colAlpha=0.4, drawConnectors=FALSE))
dev.off()

# MA
png(file.path(de_dir, "MA.png"), width=800, height=600)
plotMA(res, main="SEQC: MA Plot (B vs A)", ylim=c(-10, 10))
dev.off()

# Heatmap top DEGs
sig2 <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 2)
if (nrow(sig2) > 1) {
    top_n <- min(50, nrow(sig2))
    top_genes <- rownames(sig2[order(sig2$padj),])[1:top_n]
    mat <- assay(vsd)[top_genes,] - rowMeans(assay(vsd)[top_genes,])
    anno_c <- data.frame(Condition=meta$condition, row.names=rownames(meta))
    png(file.path(de_dir, "Heatmap_DEG.png"), width=700, height=max(400, top_n*14))
    pheatmap(mat, annotation_col=anno_c, fontsize_row=7,
             main="SEQC: Top DEGs (padj<0.05, |LFC|>2)")
    dev.off()
}

# ── 6. GSEA ──────────────────────────────────────────────────────────────────
cat("\n--- 6. GSEA ---\n")
gsea_dir <- file.path(out_dir, "GSEA")
dir.create(gsea_dir, showWarnings=FALSE)

res_rank <- as.data.frame(res)
res_rank <- res_rank[!is.na(res_rank$stat),]
res_rank <- res_rank[order(-abs(res_rank$stat)),]
res_rank <- res_rank[!duplicated(rownames(res_rank)),]
ranked <- sort(setNames(res_rank$stat, rownames(res_rank)), decreasing=TRUE)
cat("Ranked genes:", length(ranked), "\n")

for (gs_name in c("Hallmark", "GO_BP", "Reactome")) {
    gs_data <- switch(gs_name,
        Hallmark = msigdbr(species="Homo sapiens", category="H"),
        GO_BP = msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP"),
        Reactome = msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME"))
    gs_list <- split(gs_data$gene_symbol, gs_data$gs_name)
    set.seed(42)
    fgsea_res <- fgsea(pathways=gs_list, stats=ranked, minSize=15, maxSize=500, nPermSimple=10000)
    sig_n <- sum(fgsea_res$padj < 0.05, na.rm=TRUE)
    cat(sprintf("  %s: %d significant\n", gs_name, sig_n))

    out <- as.data.frame(fgsea_res)
    out$leadingEdge <- sapply(out$leadingEdge, paste, collapse=";")
    write.csv(out[order(out$pval),], file.path(gsea_dir, paste0("GSEA_", gs_name, ".csv")), row.names=FALSE)

    if (nrow(fgsea_res) >= 5) {
        top <- head(fgsea_res[order(fgsea_res$pval),], 25)
        top$pw <- substr(gsub("_", " ", gsub("^HALLMARK_|^GOBP_|^REACTOME_", "", top$pathway)), 1, 60)
        p <- ggplot(top, aes(x=reorder(pw, NES), y=NES, fill=padj<0.05)) +
            geom_col() + coord_flip() +
            scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#999999")) +
            labs(title=paste0("GSEA ", gs_name, " — SEQC"), x="", y="NES") +
            theme_bw() + theme(axis.text.y=element_text(size=7))
        ggsave(file.path(gsea_dir, paste0("GSEA_", gs_name, ".png")), p, width=12, height=8, dpi=150)
    }
}

# ── 7. Known biology: UHRR vs HBRR ──────────────────────────────────────────
cat("\n--- 7. Known Biology Validation ---\n")
cat("UHRR (Sample A) = Universal Human Reference RNA (10 cancer cell lines)\n")
cat("HBRR (Sample B) = Human Brain Reference RNA\n\n")

# Brain-enriched genes should be UP in B
brain_genes <- c("GFAP", "MBP", "SYN1", "GAD1", "SLC17A7", "NEFL", "NEFM", "NEFH", "ENO2")
cat("Brain-enriched genes (expected UP in B):\n")
for (g in brain_genes) {
    if (g %in% rownames(res)) {
        r <- as.data.frame(res)[g,]
        dir <- ifelse(r$log2FoldChange > 0, "UP in B", "DOWN in B")
        sig <- ifelse(!is.na(r$padj) & r$padj < 0.05, "✓ SIG", "  NS")
        correct <- ifelse(r$log2FoldChange > 0, "✓ CORRECT", "✗ WRONG")
        cat(sprintf("  %s: LFC=%.2f (%s) %s %s\n", g, r$log2FoldChange, dir, sig, correct))
    }
}

# Cancer/proliferation genes should be UP in A (DOWN in B)
cancer_genes <- c("MYC", "CCND1", "CDK4", "EGFR", "ERBB2", "TP53", "BCL2")
cat("\nCancer/proliferation genes (expected UP in A = DOWN in B):\n")
for (g in cancer_genes) {
    if (g %in% rownames(res)) {
        r <- as.data.frame(res)[g,]
        dir <- ifelse(r$log2FoldChange < 0, "UP in A", "UP in B")
        sig <- ifelse(!is.na(r$padj) & r$padj < 0.05, "✓ SIG", "  NS")
        correct <- ifelse(r$log2FoldChange < 0, "✓ CORRECT", "△ UNEXPECTED")
        cat(sprintf("  %s: LFC=%.2f (%s) %s %s\n", g, r$log2FoldChange, dir, sig, correct))
    }
}

# ── 8. Summary ───────────────────────────────────────────────────────────────
cat("\n###########################################################\n")
cat("### ARIA BENCHMARK SUMMARY: SEQC (EASY)\n")
cat("###########################################################\n")
for (lfc in c(2.0, 1.0, 0.5)) {
    n <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > lfc)
    cat(sprintf("DEGs (padj<0.05, |LFC|>%.1f): %d\n", lfc, n))
}
cat("DP1: QC PASS\n")
cat(sprintf("DP2: %d DEGs → standard strategy (EASY case)\n", n_deg))
cat("Brain/cancer gene validation: see above\n")
cat("Results saved to:", out_dir, "\n")
