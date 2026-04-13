"""
ssGSEA Module

Single-sample GSEA for per-sample pathway activity scoring.
"""

from pathlib import Path


class ssGSEAModule:
    """Single-sample GSEA using rank-based enrichment."""

    R_TEMPLATE = '''#!/usr/bin/env Rscript
library(DESeq2)
library(msigdbr)
library(pheatmap)
library(ggplot2)
library(reshape2)

out_dir <- "{output_dir}"
counts_file <- "{counts_file}"

# Load and normalize
counts_raw <- read.delim(counts_file, row.names=1, check.names=FALSE)
gene_names <- counts_raw$gene_name
counts_raw$gene_name <- NULL
counts_int <- round(as.matrix(counts_raw))

sample_names <- colnames(counts_int)
meta <- data.frame(
    condition = factor(ifelse(grepl("^{control_prefix}_", sample_names), "{control_label}", "{treatment_label}"),
                       levels=c("{control_label}","{treatment_label}")),
    row.names = sample_names
)

dds <- DESeqDataSetFromMatrix(counts_int, meta, ~ condition)
keep <- rowSums(counts(dds) >= 10) >= ncol(dds)/2
dds <- dds[keep,]
vsd <- vst(dds, blind=TRUE)
expr <- assay(vsd)
rownames(expr) <- gene_names[match(rownames(expr), rownames(counts_int))]
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "",]
expr <- expr[!duplicated(rownames(expr)),]

# Hallmark gene sets
hallmark <- msigdbr(species="{species}", category="H")
gs_list <- split(hallmark$gene_symbol, hallmark$gs_name)

# ssGSEA: rank-based enrichment per sample
ssgsea <- function(expr, gene_sets, min_size=10) {{
    scores <- matrix(NA, nrow=length(gene_sets), ncol=ncol(expr))
    rownames(scores) <- names(gene_sets)
    colnames(scores) <- colnames(expr)
    for (j in seq_len(ncol(expr))) {{
        ranks <- rank(expr[,j]); names(ranks) <- rownames(expr)
        for (i in seq_along(gene_sets)) {{
            gs <- intersect(gene_sets[[i]], rownames(expr))
            if (length(gs) < min_size) next
            scores[i,j] <- (mean(ranks[gs]) - mean(ranks[setdiff(names(ranks), gs)])) / length(ranks) * 2
        }}
    }}
    scores[!is.na(scores[,1]),,drop=FALSE]
}}

cat("Running ssGSEA...\\n")
ss <- ssgsea(expr, gs_list)
rownames(ss) <- gsub("_", " ", gsub("HALLMARK_", "", rownames(ss)))
write.csv(as.data.frame(ss), file.path(out_dir, "ssGSEA_scores.csv"))

# Heatmap
anno_c <- data.frame(Condition=meta$condition, row.names=rownames(meta))
anno_colors <- list(Condition=c("{control_label}"="#377EB8", "{treatment_label}"="#E41A1C"))
png(file.path(out_dir, "ssGSEA_heatmap.png"), width=900, height=900)
pheatmap(ss, annotation_col=anno_c, annotation_colors=anno_colors,
         fontsize_row=7, main="ssGSEA Hallmark Pathway Activity",
         color=colorRampPalette(c("#2166AC","white","#B2182B"))(100))
dev.off()

# t-test per pathway
tests <- data.frame()
for (i in seq_len(nrow(ss))) {{
    ctrl <- ss[i, meta$condition=="{control_label}"]
    trt <- ss[i, meta$condition=="{treatment_label}"]
    tt <- tryCatch(t.test(trt, ctrl), error=function(e) NULL)
    if (!is.null(tt))
        tests <- rbind(tests, data.frame(pathway=rownames(ss)[i],
            mean_ctrl=round(mean(ctrl),4), mean_trt=round(mean(trt),4),
            diff=round(mean(trt)-mean(ctrl),4), pvalue=signif(tt$p.value,3)))
}}
tests <- tests[order(tests$pvalue),]
write.csv(tests, file.path(out_dir, "ssGSEA_tests.csv"), row.names=FALSE)
cat("ssGSEA pathways p<0.05:", sum(tests$pvalue<0.05), "\\n")

# Key pathway boxplot
top_pw <- head(tests, 6)$pathway
if (length(top_pw) > 0) {{
    plot_data <- list()
    for (pw in top_pw) {{
        plot_data[[pw]] <- data.frame(sample=colnames(ss), score=ss[pw,],
            condition=meta$condition, pathway=pw)
    }}
    pdf <- do.call(rbind, plot_data)
    p <- ggplot(pdf, aes(x=condition, y=score, fill=condition)) +
        geom_boxplot(alpha=0.7) + geom_jitter(width=0.15, size=2) +
        facet_wrap(~pathway, scales="free_y", ncol=3) +
        scale_fill_manual(values=c("{control_label}"="#377EB8","{treatment_label}"="#E41A1C")) +
        labs(title="ssGSEA: Top Variable Pathways", y="Enrichment Score") +
        theme_bw()
    ggsave(file.path(out_dir, "ssGSEA_boxplot.png"), p, width=12, height=8, dpi=150)
}}

cat("ssGSEA complete.\\n")
'''

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_script(self, counts_file: str, species: str = "Mus musculus",
                        control_prefix: str = "WT", control_label: str = "WT",
                        treatment_label: str = "KO") -> str:
        script = self.R_TEMPLATE.format(
            output_dir=str(self.output_dir), counts_file=counts_file,
            species=species, control_prefix=control_prefix,
            control_label=control_label, treatment_label=treatment_label)
        path = self.output_dir / "run_ssGSEA.R"
        path.write_text(script)
        return str(path)
