#!/usr/bin/env Rscript
###############################################################################
# GSEA-PRIORITY BENCHMARK: Create a DEG<50 scenario
# Use Pasilla with only 2 vs 2 samples (reduced power) to simulate
# a dataset where DEGs are insufficient and GSEA becomes essential
###############################################################################

library(DESeq2)
library(fgsea)
library(msigdbr)
library(ggplot2)

out_dir <- "/data/benchmarks/results/gsea_priority"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("###########################################################\n")
cat("### GSEA-PRIORITY BENCHMARK: Low-power scenario\n")
cat("###########################################################\n")

# Load Pasilla with reduced sample size (2 vs 2)
counts <- as.matrix(read.csv("/data/benchmarks/data/pasilla/pasilla_gene_counts.tsv",
                              sep="\t", row.names="gene_id"))
meta_full <- read.csv("/data/benchmarks/data/pasilla/pasilla_sample_annotation.csv", row.names=1)
rownames(meta_full) <- sub("fb$", "", rownames(meta_full))

# Take only 2 untreated + 2 treated (simulate underpowered study)
untrt_samples <- rownames(meta_full)[meta_full$condition == "untreated"][1:2]
trt_samples <- rownames(meta_full)[meta_full$condition == "treated"][1:2]
subset_samples <- c(untrt_samples, trt_samples)

counts_sub <- counts[, subset_samples]
meta_sub <- data.frame(
    condition = factor(c("untreated","untreated","treated","treated"), levels=c("untreated","treated")),
    row.names = subset_samples
)

cat("Samples:", paste(subset_samples, collapse=", "), "\n")
cat("Design: 2 untreated vs 2 treated (deliberately underpowered)\n\n")

# ── DE Analysis ──
keep <- rowSums(counts_sub >= 10) >= 2
counts_filt <- counts_sub[keep, ]
cat("Genes after filtering:", nrow(counts_filt), "\n")

dds <- DESeqDataSetFromMatrix(counts_filt, meta_sub, ~ condition)
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds, contrast=c("condition","treated","untreated"), alpha=0.05)

for (lfc in c(1.0, 0.5, 0.0)) {
    n <- sum(!is.na(res$padj) & res$padj<0.05 & abs(res$log2FoldChange)>lfc)
    cat(sprintf("  padj<0.05, |LFC|>%.1f: %d DEGs\n", lfc, n))
}

n_deg1 <- sum(!is.na(res$padj) & res$padj<0.05 & abs(res$log2FoldChange)>1)

# ── DP2 Decision ──
cat("\n--- DP2 Decision ---\n")
if (n_deg1 >= 100) {
    strategy <- "standard"
} else if (n_deg1 >= 10) {
    strategy <- "augmented"
} else {
    strategy <- "gsea_priority"
}
cat(sprintf("DEGs (|LFC|>1) = %d → Strategy: %s\n", n_deg1, strategy))

if (strategy == "gsea_priority") {
    cat("\n*** GSEA-PRIORITY MODE TRIGGERED ***\n")
    cat("With <10 DEGs, ORA is underpowered.\n")
    cat("GSEA uses full gene ranking — does not depend on DEG cutoff.\n\n")
}

# ── GSEA (the key analysis in gsea_priority mode) ──
cat("--- GSEA (critical in low-DEG scenario) ---\n")
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$stat),]
ranked <- sort(setNames(res_df$stat, rownames(res_df)), decreasing=TRUE)
cat("Ranked genes:", length(ranked), "\n")

# Load Drosophila gene sets (Pasilla is Drosophila)
# MSigDB doesn't have Drosophila directly, so use orthologs or GO terms
# Use GO terms mapped to Drosophila gene IDs
# Actually, fgsea works with any gene ID as long as they match
# For this test, we'll use the gene IDs directly with generic pathway sets

# Create simple gene sets from the data itself (top/bottom ranked gene groups)
# This tests whether GSEA can detect directional trends even with n=2

# Better approach: use msigdbr with Drosophila
dmel_go <- tryCatch(
    msigdbr(species="Drosophila melanogaster", category="C5", subcategory="GO:BP"),
    error = function(e) NULL
)

if (!is.null(dmel_go) && nrow(dmel_go) > 0) {
    gs_list <- split(dmel_go$gene_symbol, dmel_go$gs_name)
    cat("Drosophila GO:BP gene sets:", length(gs_list), "\n")
} else {
    cat("Drosophila gene sets not available in msigdbr, using manual approach\n")
    # Create gene sets from FlyBase GO annotations
    # For this benchmark, we'll demonstrate the DP2 logic works
    gs_list <- list()
}

if (length(gs_list) > 0) {
    set.seed(42)
    gsea_res <- fgsea(pathways=gs_list, stats=ranked, minSize=10, maxSize=500, nPermSimple=10000)
    gsea_sig <- sum(gsea_res$padj < 0.05, na.rm=TRUE)
    cat(sprintf("GSEA GO:BP significant: %d pathways\n", gsea_sig))

    gsea_out <- as.data.frame(gsea_res)
    gsea_out$leadingEdge <- sapply(gsea_out$leadingEdge, paste, collapse=";")
    write.csv(gsea_out[order(gsea_out$pval),], file.path(out_dir, "GSEA_GO_BP.csv"), row.names=FALSE)
} else {
    gsea_sig <- NA
    cat("GSEA skipped (no gene sets available for Drosophila)\n")
}

# ── Comparison: Full power (4 vs 3) vs Reduced power (2 vs 2) ──
cat("\n--- Power comparison ---\n")
# Full power results (from previous Pasilla benchmark)
cat("Full power (4 vs 3): 224 DEGs (|LFC|>1), 635 DEGs (|LFC|>0.5)\n")
cat(sprintf("Reduced power (2 vs 2): %d DEGs (|LFC|>1)\n", n_deg1))
cat(sprintf("DEG loss: %.0f%%\n", (1 - n_deg1/224) * 100))
cat(sprintf("DP2 strategy change: 'standard' → '%s'\n", strategy))

# ── Summary ──
cat("\n###########################################################\n")
cat("### GSEA-PRIORITY SUMMARY\n")
cat("###########################################################\n")
summary_df <- data.frame(
    Metric = c("Samples", "DEGs (|LFC|>1)", "DEGs (|LFC|>0.5)", "DP2 Strategy",
               "GSEA pathways (padj<0.05)", "ORA feasible?"),
    Full_Power = c("4 vs 3", "224", "635", "standard", "N/A (not tested)", "Yes"),
    Reduced_Power = c("2 vs 2", as.character(n_deg1),
                      as.character(sum(!is.na(res$padj) & res$padj<0.05 & abs(res$log2FoldChange)>0.5)),
                      strategy,
                      ifelse(is.na(gsea_sig), "N/A", as.character(gsea_sig)),
                      ifelse(n_deg1 >= 10, "Marginal", "No — too few DEGs"))
)
print(summary_df)
write.csv(summary_df, file.path(out_dir, "gsea_priority_summary.csv"), row.names=FALSE)
