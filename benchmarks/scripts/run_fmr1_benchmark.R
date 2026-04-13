#!/usr/bin/env Rscript
###############################################################################
# ARIA Benchmark: Fmr1 KO (GSE180135) — HARD scenario
# Mouse cortical neurons, Fmr1 KO vs WT, n=3 vs 3
# Tests: DP1 (QC), DP2 (few DEGs → GSEA), DP4 (signatures), DP5 (validation)
###############################################################################

library(DESeq2)
library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(futile.logger)

out_dir <- "/data/benchmarks/results/fmr1"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("###########################################################\n")
cat("### ARIA Benchmark: Fmr1 KO (GSE180135) — HARD\n")
cat("###########################################################\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
cat("\n--- 1. Data Loading ---\n")
counts <- read.delim("/data/benchmarks/data/fmr1/GSE180135_deSeq2_counts.txt",
                      row.names=1, check.names=FALSE)
counts <- as.matrix(counts)
cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Samples:", paste(colnames(counts), collapse=", "), "\n")

meta <- data.frame(
    genotype = factor(ifelse(grepl("^WT", colnames(counts)), "WT", "KO"), levels=c("WT","KO")),
    row.names = colnames(counts)
)
print(meta)

# ── 2. DP1: QC Assessment ───────────────────────────────────────────────────
cat("\n--- 2. DP1: QC Assessment ---\n")
total_counts <- colSums(counts)
cat("Total counts per sample:\n"); print(total_counts)
cat("Non-zero genes:", sum(rowSums(counts) > 0), "\n")
cat(sprintf("DP1 DECISION: Min counts = %d → %s\n", min(total_counts),
    ifelse(min(total_counts) > 1e6, "PASS", "WARNING")))

# ── 3. DE Analysis ──────────────────────────────────────────────────────────
cat("\n--- 3. DE Analysis (DESeq2) ---\n")
de_dir <- file.path(out_dir, "DE")
dir.create(de_dir, showWarnings=FALSE)

keep <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep, ]
cat("Genes after filtering:", nrow(counts_filt), "\n")

dds <- DESeqDataSetFromMatrix(counts_filt, meta, ~ genotype)
dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype","KO","WT"), alpha=0.05)

# Gene symbols (already gene names as rownames)
res$symbol <- rownames(res)

# Map to ENTREZ
entrez <- mapIds(org.Mm.eg.db, keys=rownames(res), column="ENTREZID",
                  keytype="SYMBOL", multiVals="first")
res$entrez <- entrez[rownames(res)]

for (lfc in c(1.0, 0.5, 0.0)) {
    n <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > lfc)
    up <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > lfc)
    down <- sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -lfc)
    cat(sprintf("  padj<0.05, |LFC|>%.1f: %d DEGs (Up: %d, Down: %d)\n", lfc, n, up, down))
}

# Save
res_df <- as.data.frame(res)
write.csv(res_df, file.path(de_dir, "DE_all.csv"), row.names=TRUE)

sig1 <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
sig05 <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5)
write.csv(sig1[order(sig1$padj),], file.path(de_dir, "DE_sig_LFC1.csv"), row.names=TRUE)
write.csv(sig05[order(sig05$padj),], file.path(de_dir, "DE_sig_LFC05.csv"), row.names=TRUE)

# ── 4. DP2: DE Evaluation — KEY DECISION POINT ─────────────────────────────
cat("\n--- 4. DP2: DE Result Evaluation (KEY) ---\n")
n_deg1 <- nrow(sig1)
n_deg05 <- nrow(sig05)

if (n_deg1 >= 100) {
    cat(sprintf("DP2 DECISION: %d DEGs (|LFC|>1) → Sufficient for ORA + GSEA\n", n_deg1))
    strategy <- "standard"
} else if (n_deg1 >= 10) {
    cat(sprintf("DP2 DECISION: %d DEGs (|LFC|>1) → ADD GSEA, TRY relaxed cutoffs\n", n_deg1))
    strategy <- "augmented"
} else {
    cat(sprintf("DP2 DECISION: %d DEGs (|LFC|>1) → FEW DEGs! PRIORITIZE GSEA over ORA\n", n_deg1))
    cat("  → GSEA uses full gene ranking, not dependent on arbitrary DEG cutoff\n")
    cat("  → This is the HARD scenario where ARIA's adaptive strategy is critical\n")
    strategy <- "gsea_priority"
}
cat("Strategy:", strategy, "\n")

# ── 5. Plots ─────────────────────────────────────────────────────────────────
cat("\n--- 5. Visualization ---\n")

vsd <- vst(dds, blind=FALSE)

# PCA
pca_data <- plotPCA(vsd, intgroup="genotype", returnData=TRUE)
pct_var <- round(100*attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color=genotype, label=name)) +
    geom_point(size=5) + geom_text_repel(size=3.5) +
    xlab(paste0("PC1: ", pct_var[1], "%")) + ylab(paste0("PC2: ", pct_var[2], "%")) +
    ggtitle("Fmr1 KO: PCA (KO vs WT cortical neurons)") + theme_bw() +
    scale_color_manual(values=c("WT"="#377EB8","KO"="#E41A1C"))
ggsave(file.path(de_dir, "PCA.png"), p, width=8, height=6, dpi=150)

# Volcano
for (lfc_cut in c(1.0, 0.5)) {
    lfc_label <- gsub("\\.", "", as.character(lfc_cut))
    png(file.path(de_dir, paste0("Volcano_LFC", lfc_label, ".png")), width=900, height=700)
    print(EnhancedVolcano(res, lab=res$symbol, x="log2FoldChange", y="padj",
        pCutoff=0.05, FCcutoff=lfc_cut,
        title=paste0("Fmr1 KO vs WT (|LFC|>", lfc_cut, ")"),
        pointSize=2, labSize=3.5, colAlpha=0.6))
    dev.off()
}

# Heatmap
if (nrow(sig05) > 1) {
    top_n <- min(50, nrow(sig05))
    top_genes <- rownames(sig05[order(sig05$padj),])[1:top_n]
    mat <- assay(vsd)[top_genes,] - rowMeans(assay(vsd)[top_genes,])
    anno_c <- data.frame(Genotype=meta$genotype, row.names=rownames(meta))
    png(file.path(de_dir, "Heatmap_DEG.png"), width=700, height=max(400, top_n*14))
    pheatmap(mat, annotation_col=anno_c, fontsize_row=7,
             main="Fmr1 KO: Top DEGs (padj<0.05, |LFC|>0.5)")
    dev.off()
}

# ── 6. DP5: Cross-method validation ─────────────────────────────────────────
cat("\n--- 6. DP5: Cross-Method Validation ---\n")
dp5_dir <- file.path(out_dir, "method_comparison")
dir.create(dp5_dir, showWarnings=FALSE)

y <- DGEList(counts=counts_filt, group=meta$genotype)
y <- y[filterByExpr(y, group=meta$genotype),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ meta$genotype)
y <- estimateDisp(y, design)

# edgeR exact
et <- exactTest(y, pair=c("WT","KO"))
res_exact <- topTags(et, n=Inf)$table

# limma-voom
v <- voom(y, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef=2, number=Inf)

cat("Cross-method (padj<0.05, |LFC|>0.5):\n")
n_deseq2 <- nrow(sig05)
n_edger <- sum(!is.na(res_exact$FDR) & res_exact$FDR < 0.05 & abs(res_exact$logFC) > 0.5)
n_limma <- sum(!is.na(res_limma$adj.P.Val) & res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > 0.5)
cat(sprintf("  DESeq2: %d, edgeR exact: %d, limma-voom: %d\n", n_deseq2, n_edger, n_limma))

# LFC correlation
common <- intersect(rownames(res), rownames(res_exact))
cor_val <- cor(as.data.frame(res)[common,"log2FoldChange"], res_exact[common,"logFC"], use="complete.obs")
cat(sprintf("  LFC correlation (DESeq2 vs edgeR): r=%.4f\n", cor_val))

# ── 7. GSEA (critical for HARD scenario) ────────────────────────────────────
cat("\n--- 7. GSEA (CRITICAL for hard scenario) ---\n")
gsea_dir <- file.path(out_dir, "GSEA")
dir.create(gsea_dir, showWarnings=FALSE)

res_rank <- as.data.frame(res)
res_rank <- res_rank[!is.na(res_rank$stat),]
res_rank <- res_rank[order(-abs(res_rank$stat)),]
res_rank <- res_rank[!duplicated(rownames(res_rank)),]
ranked <- sort(setNames(res_rank$stat, rownames(res_rank)), decreasing=TRUE)
cat("Ranked genes:", length(ranked), "\n")

collections <- list(
    Hallmark = msigdbr(species="Mus musculus", category="H"),
    GO_BP = msigdbr(species="Mus musculus", category="C5", subcategory="GO:BP"),
    GO_CC = msigdbr(species="Mus musculus", category="C5", subcategory="GO:CC"),
    KEGG = msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_MEDICUS"),
    Reactome = msigdbr(species="Mus musculus", category="C2", subcategory="CP:REACTOME")
)

gsea_summary <- data.frame()
for (gs_name in names(collections)) {
    gs_list <- split(collections[[gs_name]]$gene_symbol, collections[[gs_name]]$gs_name)
    set.seed(42)
    fgsea_res <- fgsea(pathways=gs_list, stats=ranked, minSize=15, maxSize=500, nPermSimple=10000)
    sig_n <- sum(fgsea_res$padj < 0.05, na.rm=TRUE)
    cat(sprintf("  %s: %d significant (padj<0.05)\n", gs_name, sig_n))
    gsea_summary <- rbind(gsea_summary, data.frame(Collection=gs_name, Significant=sig_n))

    out <- as.data.frame(fgsea_res)
    out$leadingEdge <- sapply(out$leadingEdge, paste, collapse=";")
    write.csv(out[order(out$pval),], file.path(gsea_dir, paste0("GSEA_", gs_name, ".csv")), row.names=FALSE)

    if (nrow(fgsea_res) >= 5) {
        top <- head(fgsea_res[order(fgsea_res$pval),], 25)
        top$pw <- substr(gsub("_", " ", gsub("^HALLMARK_|^GOBP_|^GOCC_|^KEGG_MEDICUS_|^REACTOME_", "", top$pathway)), 1, 60)
        p <- ggplot(top, aes(x=reorder(pw, NES), y=NES, fill=padj<0.05)) +
            geom_col() + coord_flip() +
            scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#999999")) +
            labs(title=paste0("GSEA ", gs_name, " — Fmr1 KO"), x="", y="NES") +
            theme_bw() + theme(axis.text.y=element_text(size=7))
        ggsave(file.path(gsea_dir, paste0("GSEA_", gs_name, ".png")), p, width=12, height=8, dpi=150)
    }
}

# ── 8. Known Fmr1 biology validation ────────────────────────────────────────
cat("\n--- 8. Known Fmr1/FMRP Target Validation ---\n")
cat("Known FMRP-regulated genes / Fragile X biology:\n")
known_genes <- c("Fmr1", "Map1b", "Psd95", "Dlg4", "Shank3", "Arc", "Camk2a",
                  "Mmp9", "Nlgn1", "Grin1", "Grin2a", "Grin2b", "Bdnf")
for (g in known_genes) {
    if (g %in% rownames(res)) {
        r <- as.data.frame(res)[g,]
        status <- ifelse(!is.na(r$padj) & r$padj < 0.05, "✓ SIG", "  NS")
        cat(sprintf("  %s: LFC=%.3f, padj=%s %s\n", g, r$log2FoldChange,
            ifelse(is.na(r$padj), "NA", sprintf("%.1e", r$padj)), status))
    } else {
        cat(sprintf("  %s: NOT FOUND\n", g))
    }
}

# ── 9. Summary ───────────────────────────────────────────────────────────────
cat("\n###########################################################\n")
cat("### ARIA BENCHMARK SUMMARY: Fmr1 KO (HARD)\n")
cat("###########################################################\n")
cat(sprintf("Genes tested: %d\n", nrow(res)))
cat(sprintf("DEGs (padj<0.05, |LFC|>1): %d\n", nrow(sig1)))
cat(sprintf("DEGs (padj<0.05, |LFC|>0.5): %d\n", nrow(sig05)))
cat(sprintf("DEGs (padj<0.05, no LFC): %d\n", sum(!is.na(res$padj) & res$padj<0.05)))
cat(sprintf("Strategy selected: %s\n", strategy))
cat("\nGSEA results:\n")
print(gsea_summary)
cat(sprintf("\nDP1: QC %s\n", ifelse(min(total_counts)>1e6, "PASS", "WARNING")))
cat(sprintf("DP2: %d DEGs → %s strategy\n", n_deg1, strategy))
cat(sprintf("DP5: Cross-validation LFC r=%.4f\n", cor_val))
cat("Results saved to:", out_dir, "\n")
