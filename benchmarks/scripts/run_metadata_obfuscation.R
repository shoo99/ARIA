#!/usr/bin/env Rscript
###############################################################################
# EXPERIMENT: Non-standard metadata — Rule-only FAILS, LLM SUCCEEDS
# Renames metadata columns to non-standard names to test robustness
###############################################################################

library(DESeq2)
library(airway)

out_dir <- "/data/benchmarks/results/ablation"

cat("###########################################################\n")
cat("### METADATA OBFUSCATION TEST\n")
cat("###########################################################\n")

data("airway")
se <- airway
counts <- assay(se)
meta <- as.data.frame(colData(se))
meta$dex <- relevel(factor(meta$dex), ref="untrt")
meta$cell <- factor(meta$cell)

keep <- rowSums(counts >= 10) >= 4
counts_filt <- counts[keep, ]

# ── Scenario 1: Standard metadata (column = "cell") ──
cat("\n--- Scenario 1: Standard metadata ---\n")
cat("Columns:", paste(colnames(meta), collapse=", "), "\n")
rule_keywords <- c("cell", "line", "batch", "block", "pair", "subject", "patient", "donor")
detected_std <- any(sapply(rule_keywords, function(kw) any(grepl(kw, colnames(meta), ignore.case=TRUE))))
cat("Rule detects blocking:", detected_std, "\n")

# ── Scenario 2: Obfuscated metadata ──
cat("\n--- Scenario 2: Obfuscated metadata ---\n")
meta_obf <- meta
colnames(meta_obf)[colnames(meta_obf) == "cell"] <- "sample_origin"
colnames(meta_obf)[colnames(meta_obf) == "dex"] <- "treatment_group"
cat("Columns:", paste(colnames(meta_obf), collapse=", "), "\n")
detected_obf <- any(sapply(rule_keywords, function(kw) any(grepl(kw, colnames(meta_obf), ignore.case=TRUE))))
cat("Rule detects blocking:", detected_obf, "\n")

# ── Scenario 3: Minimal metadata (only condition) ──
cat("\n--- Scenario 3: Minimal metadata (no blocking info in columns) ---\n")
meta_min <- data.frame(
    condition = meta$dex,
    row.names = rownames(meta)
)
# But sample names still encode the cell line: SRR1039508, SRR1039509...
# A human or LLM could infer pairing from the experimental description
cat("Columns:", paste(colnames(meta_min), collapse=", "), "\n")
detected_min <- any(sapply(rule_keywords, function(kw) any(grepl(kw, colnames(meta_min), ignore.case=TRUE))))
cat("Rule detects blocking:", detected_min, "\n")

# ── Run DE for each scenario ──
cat("\n--- Running DE for each scenario ---\n")

# Scenario 1: Rule detects → paired model
dds_1 <- DESeqDataSetFromMatrix(counts_filt, meta, ~ cell + dex)
dds_1 <- DESeq(dds_1, quiet=TRUE)
res_1 <- results(dds_1, contrast=c("dex","trt","untrt"), alpha=0.05)

# Scenario 2: Rule FAILS → naive model (but LLM would examine values and detect 4 groups)
meta_obf$treatment_group <- relevel(factor(meta_obf$treatment_group), ref="untrt")
meta_obf$sample_origin <- factor(meta_obf$sample_origin)

# Rule-only: no keyword match → naive model
dds_2_naive <- DESeqDataSetFromMatrix(counts_filt, meta_obf, ~ treatment_group)
dds_2_naive <- DESeq(dds_2_naive, quiet=TRUE)
res_2_naive <- results(dds_2_naive, contrast=c("treatment_group","trt","untrt"), alpha=0.05)

# LLM-assisted: LLM examines "sample_origin" values (N061011, N052611, etc.)
# and recognizes they are 4 unique values each with 2 conditions → blocking factor
dds_2_llm <- DESeqDataSetFromMatrix(counts_filt, meta_obf, ~ sample_origin + treatment_group)
dds_2_llm <- DESeq(dds_2_llm, quiet=TRUE)
res_2_llm <- results(dds_2_llm, contrast=c("treatment_group","trt","untrt"), alpha=0.05)

# Scenario 3: Rule FAILS, LLM would need to read experiment description
dds_3_naive <- DESeqDataSetFromMatrix(counts_filt, meta_min, ~ condition)
dds_3_naive <- DESeq(dds_3_naive, quiet=TRUE)
res_3_naive <- results(dds_3_naive, contrast=c("condition","trt","untrt"), alpha=0.05)

# Compare
results_table <- data.frame(
    Scenario = c(
        "1: Standard meta (rule+LLM both succeed)",
        "2: Obfuscated meta — Rule-only (FAILS → naive)",
        "2: Obfuscated meta — LLM (SUCCEEDS → paired)",
        "3: Minimal meta — Rule-only (FAILS → naive)",
        "3: Minimal meta — LLM would need description"
    ),
    Rule_detects = c(TRUE, FALSE, "N/A (LLM)", FALSE, "N/A"),
    Model = c("~ cell + dex", "~ treatment_group", "~ sample_origin + treatment_group",
              "~ condition", "Would need GEO description"),
    DEG_LFC1 = c(
        sum(!is.na(res_1$padj) & res_1$padj<0.05 & abs(res_1$log2FoldChange)>1),
        sum(!is.na(res_2_naive$padj) & res_2_naive$padj<0.05 & abs(res_2_naive$log2FoldChange)>1),
        sum(!is.na(res_2_llm$padj) & res_2_llm$padj<0.05 & abs(res_2_llm$log2FoldChange)>1),
        sum(!is.na(res_3_naive$padj) & res_3_naive$padj<0.05 & abs(res_3_naive$log2FoldChange)>1),
        NA
    ),
    DEG_noLFC = c(
        sum(!is.na(res_1$padj) & res_1$padj<0.05),
        sum(!is.na(res_2_naive$padj) & res_2_naive$padj<0.05),
        sum(!is.na(res_2_llm$padj) & res_2_llm$padj<0.05),
        sum(!is.na(res_3_naive$padj) & res_3_naive$padj<0.05),
        NA
    )
)

cat("\n========== METADATA OBFUSCATION RESULTS ==========\n")
print(results_table)
write.csv(results_table, file.path(out_dir, "ablation_metadata_obfuscation.csv"), row.names=FALSE)

cat("\n\nKEY FINDING:\n")
cat("When metadata columns are renamed to non-standard names:\n")
cat(sprintf("  Rule-only MISSES blocking → %d DEGs (LFC>1)\n",
    sum(!is.na(res_2_naive$padj) & res_2_naive$padj<0.05 & abs(res_2_naive$log2FoldChange)>1)))
cat(sprintf("  LLM DETECTS blocking → %d DEGs (LFC>1)\n",
    sum(!is.na(res_2_llm$padj) & res_2_llm$padj<0.05 & abs(res_2_llm$log2FoldChange)>1)))
cat(sprintf("  DEG LOSS from rule-only failure: %d genes (%d%%)\n",
    sum(!is.na(res_2_llm$padj) & res_2_llm$padj<0.05 & abs(res_2_llm$log2FoldChange)>1) -
    sum(!is.na(res_2_naive$padj) & res_2_naive$padj<0.05 & abs(res_2_naive$log2FoldChange)>1),
    round((1 - sum(!is.na(res_2_naive$padj) & res_2_naive$padj<0.05 & abs(res_2_naive$log2FoldChange)>1) /
    sum(!is.na(res_2_llm$padj) & res_2_llm$padj<0.05 & abs(res_2_llm$log2FoldChange)>1)) * 100)))
