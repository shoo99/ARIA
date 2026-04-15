#!/usr/bin/env Rscript
###############################################################################
# REPRODUCIBILITY TEST: 3x independent runs on Airway
# Tests whether statistical results are identical across runs
###############################################################################

library(DESeq2)
library(fgsea)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(airway)

out_dir <- "/data/benchmarks/results/reproducibility"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("###########################################################\n")
cat("### REPRODUCIBILITY TEST: 3x Airway Benchmark\n")
cat("###########################################################\n")

data("airway")
se <- airway
counts <- assay(se)
meta <- as.data.frame(colData(se))
meta$dex <- relevel(factor(meta$dex), ref="untrt")
meta$cell <- factor(meta$cell)

keep <- rowSums(counts >= 10) >= 4
counts_filt <- counts[keep, ]

# Hallmark gene sets
hallmark <- msigdbr(species="Homo sapiens", category="H")
gs_list <- split(hallmark$gene_symbol, hallmark$gs_name)

###############################################################################
# Run 3 independent analyses
###############################################################################

results_list <- list()
gsea_list <- list()
deg_counts <- data.frame()

for (run in 1:3) {
    cat(sprintf("\n--- Run %d ---\n", run))

    # DE analysis (deterministic with same seed)
    set.seed(42)  # Same seed each run
    dds <- DESeqDataSetFromMatrix(counts_filt, meta, ~ cell + dex)
    dds <- DESeq(dds, quiet=TRUE)
    res <- results(dds, contrast=c("dex","trt","untrt"), alpha=0.05)

    n1 <- sum(!is.na(res$padj) & res$padj<0.05 & abs(res$log2FoldChange)>1)
    n05 <- sum(!is.na(res$padj) & res$padj<0.05 & abs(res$log2FoldChange)>0.5)
    n0 <- sum(!is.na(res$padj) & res$padj<0.05)

    cat(sprintf("  DEGs: LFC>1=%d, LFC>0.5=%d, noLFC=%d\n", n1, n05, n0))
    deg_counts <- rbind(deg_counts, data.frame(Run=run, LFC1=n1, LFC05=n05, noLFC=n0))
    results_list[[run]] <- as.data.frame(res)

    # GSEA (with fixed seed — should be deterministic)
    res_df <- as.data.frame(res)
    res_df$symbol <- mapIds(org.Hs.eg.db, keys=rownames(res_df),
                             column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    res_df <- res_df[!is.na(res_df$symbol) & !is.na(res_df$stat),]
    res_df <- res_df[order(-abs(res_df$stat)),]
    res_df <- res_df[!duplicated(res_df$symbol),]
    ranked <- sort(setNames(res_df$stat, res_df$symbol), decreasing=TRUE)

    set.seed(42)  # Same seed for GSEA
    gsea_res <- fgsea(pathways=gs_list, stats=ranked, minSize=15, maxSize=500, nPermSimple=10000)
    gsea_sig <- sum(gsea_res$padj < 0.05, na.rm=TRUE)
    cat(sprintf("  GSEA Hallmark significant: %d\n", gsea_sig))
    gsea_list[[run]] <- gsea_res
}

###############################################################################
# Compare runs
###############################################################################
cat("\n========== REPRODUCIBILITY RESULTS ==========\n")

cat("\nDEG counts across runs:\n")
print(deg_counts)

# Check if DEG counts are identical
deg_identical <- all(deg_counts$LFC1 == deg_counts$LFC1[1]) &&
                 all(deg_counts$LFC05 == deg_counts$LFC05[1]) &&
                 all(deg_counts$noLFC == deg_counts$noLFC[1])
cat(sprintf("\nDEG counts identical across runs: %s\n", deg_identical))

# Compare LFC values
lfc_cors <- c()
for (i in 1:2) {
    for (j in (i+1):3) {
        common <- intersect(rownames(results_list[[i]]), rownames(results_list[[j]]))
        r <- cor(results_list[[i]][common, "log2FoldChange"],
                 results_list[[j]][common, "log2FoldChange"], use="complete.obs")
        lfc_cors <- c(lfc_cors, r)
        cat(sprintf("LFC correlation Run %d vs Run %d: r=%.10f\n", i, j, r))
    }
}

# Compare padj values
padj_cors <- c()
for (i in 1:2) {
    for (j in (i+1):3) {
        common <- intersect(rownames(results_list[[i]]), rownames(results_list[[j]]))
        v1 <- results_list[[i]][common, "padj"]
        v2 <- results_list[[j]][common, "padj"]
        valid <- !is.na(v1) & !is.na(v2)
        r <- cor(v1[valid], v2[valid])
        padj_cors <- c(padj_cors, r)
        cat(sprintf("padj correlation Run %d vs Run %d: r=%.10f\n", i, j, r))
    }
}

# Compare GSEA NES
gsea_cors <- c()
for (i in 1:2) {
    for (j in (i+1):3) {
        common_pw <- intersect(gsea_list[[i]]$pathway, gsea_list[[j]]$pathway)
        nes_i <- gsea_list[[i]]$NES[match(common_pw, gsea_list[[i]]$pathway)]
        nes_j <- gsea_list[[j]]$NES[match(common_pw, gsea_list[[j]]$pathway)]
        r <- cor(nes_i, nes_j, use="complete.obs")
        gsea_cors <- c(gsea_cors, r)
        cat(sprintf("GSEA NES correlation Run %d vs Run %d: r=%.6f\n", i, j, r))
    }
}

# GSEA significant pathway overlap
gsea_sigs <- lapply(gsea_list, function(g) g$pathway[g$padj < 0.05])
overlap_12 <- length(intersect(gsea_sigs[[1]], gsea_sigs[[2]])) / length(union(gsea_sigs[[1]], gsea_sigs[[2]]))
overlap_13 <- length(intersect(gsea_sigs[[1]], gsea_sigs[[3]])) / length(union(gsea_sigs[[1]], gsea_sigs[[3]]))
overlap_23 <- length(intersect(gsea_sigs[[2]], gsea_sigs[[3]])) / length(union(gsea_sigs[[2]], gsea_sigs[[3]]))
cat(sprintf("\nGSEA significant pathway Jaccard overlap:\n"))
cat(sprintf("  Run 1 vs 2: %.4f\n", overlap_12))
cat(sprintf("  Run 1 vs 3: %.4f\n", overlap_13))
cat(sprintf("  Run 2 vs 3: %.4f\n", overlap_23))

###############################################################################
# Summary
###############################################################################
summary_df <- data.frame(
    Metric = c("DEG counts (LFC>1)", "DEG counts (LFC>0.5)", "DEG counts (no LFC)",
               "LFC correlation (min)", "padj correlation (min)",
               "GSEA NES correlation (min)", "GSEA pathway Jaccard (min)"),
    Run1_vs_Run2_vs_Run3 = c(
        paste(deg_counts$LFC1, collapse="/"),
        paste(deg_counts$LFC05, collapse="/"),
        paste(deg_counts$noLFC, collapse="/"),
        sprintf("%.10f", min(lfc_cors)),
        sprintf("%.10f", min(padj_cors)),
        sprintf("%.6f", min(gsea_cors)),
        sprintf("%.4f", min(overlap_12, overlap_13, overlap_23))
    ),
    Identical = c(
        deg_identical, deg_identical, deg_identical,
        min(lfc_cors) == 1.0, min(padj_cors) > 0.9999,
        min(gsea_cors) > 0.99, min(overlap_12, overlap_13, overlap_23) > 0.9
    )
)
print(summary_df)
write.csv(summary_df, file.path(out_dir, "reproducibility_summary.csv"), row.names=FALSE)

cat("\n###########################################################\n")
cat("### CONCLUSION\n")
cat("###########################################################\n")
cat("Statistical components (DESeq2, fgsea with fixed seed):\n")
cat(sprintf("  DEG: %s identical\n", ifelse(deg_identical, "FULLY", "NOT")))
cat(sprintf("  LFC: r=%.10f\n", min(lfc_cors)))
cat(sprintf("  GSEA NES: r=%.6f\n", min(gsea_cors)))
cat("\nLLM components (DP6, DP8): NOT tested here — require API calls\n")
cat("  → Planned for separate evaluation with BLEU/cosine similarity\n")
cat("\nResults saved to:", out_dir, "\n")
