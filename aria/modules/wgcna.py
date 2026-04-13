"""
WGCNA Module

Weighted Gene Co-expression Network Analysis.
"""

from pathlib import Path


class WGCNAModule:
    """WGCNA co-expression network analysis."""

    R_TEMPLATE = '''#!/usr/bin/env Rscript
library(WGCNA)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)

options(stringsAsFactors=FALSE)
allowWGCNAThreads(4)
cor <- WGCNA::cor

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
expr <- t(assay(vsd))

# Top variable genes
gene_var <- apply(assay(vsd), 1, var)
top_idx <- order(gene_var, decreasing=TRUE)[1:min({n_genes}, length(gene_var))]
expr <- expr[, top_idx]
cat("Using", ncol(expr), "genes,", nrow(expr), "samples\\n")

# Soft threshold
powers <- c(1:20)
sft <- pickSoftThreshold(expr, powerVector=powers, verbose=3, networkType="signed")
power <- sft$powerEstimate
if (is.na(power) || power < 1) power <- 12

png(file.path(out_dir, "soft_threshold.png"), width=900, height=450)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold", ylab="Scale Free R²", type="n", main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=0.9, col="red")
abline(h=0.85, col="red", lty=2)
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold", ylab="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

cat("Power:", power, "\\n")

# Build network
net <- blockwiseModules(expr, power=power, networkType="signed", TOMType="signed",
    minModuleSize=30, reassignThreshold=0, mergeCutHeight=0.25,
    numericLabels=TRUE, pamRespectsDendro=FALSE, saveTOMs=FALSE, verbose=3)

cor <- stats::cor
module_colors <- labels2colors(net$colors)
n_modules <- length(unique(module_colors)) - ("grey" %in% module_colors)
cat("Modules:", n_modules, "\\n")

png(file.path(out_dir, "module_dendrogram.png"), width=1200, height=600)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
    "Module Colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, main="Gene Dendrogram")
dev.off()

# Module-trait correlation
traits <- data.frame(condition_trt=as.numeric(meta$condition=="{treatment_label}"), row.names=rownames(meta))
MEs <- moduleEigengenes(expr, colors=module_colors)$eigengenes
MEs <- orderMEs(MEs)

mod_cor <- cor(MEs, traits, use="p")
mod_pval <- corPvalueStudent(mod_cor, nrow(expr))

text_mat <- paste(signif(mod_cor,2), "\\n(", signif(mod_pval,1), ")", sep="")
dim(text_mat) <- dim(mod_cor)

png(file.path(out_dir, "module_trait_correlation.png"), width=500, height=max(400, ncol(MEs)*30))
par(mar=c(6,10,3,2))
labeledHeatmap(Matrix=mod_cor, xLabels=colnames(traits), yLabels=colnames(MEs),
    colorLabels=FALSE, colors=blueWhiteRed(50), textMatrix=text_mat,
    setStdMargins=FALSE, cex.text=0.7, zlim=c(-1,1), main="Module-Trait Correlations")
dev.off()

# Save results
geno_results <- data.frame(Module=gsub("ME","",names(mod_cor[,1])),
    Correlation=round(mod_cor[,1],4), Pvalue=signif(mod_pval[,1],3),
    Size=sapply(gsub("ME","",names(mod_cor[,1])), function(c) sum(module_colors==c)))
geno_results <- geno_results[order(geno_results$Pvalue),]
write.csv(geno_results, file.path(out_dir, "module_correlations.csv"), row.names=FALSE)

# Module eigengene heatmap
png(file.path(out_dir, "eigengene_heatmap.png"), width=800, height=max(400, nrow(t(MEs))*25))
anno_c <- data.frame(Condition=meta$condition, row.names=rownames(meta))
pheatmap(t(as.matrix(MEs)), annotation_col=anno_c, main="Module Eigengenes")
dev.off()

cat("WGCNA complete. Modules:", n_modules, "\\n")
print(geno_results)
'''

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_script(self, counts_file: str, n_genes: int = 5000,
                        control_prefix: str = "WT", control_label: str = "WT",
                        treatment_label: str = "KO") -> str:
        script = self.R_TEMPLATE.format(
            output_dir=str(self.output_dir), counts_file=counts_file,
            n_genes=n_genes, control_prefix=control_prefix,
            control_label=control_label, treatment_label=treatment_label)
        path = self.output_dir / "run_WGCNA.R"
        path.write_text(script)
        return str(path)
