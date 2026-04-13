"""
Cell Type Deconvolution Module

Marker-based cell type composition estimation from bulk RNA-seq.
"""

from pathlib import Path


class DeconvolutionModule:
    """Cell type deconvolution using marker genes."""

    # Default brain cell type markers
    BRAIN_MARKERS = {
        "Neurons_Excitatory": ["Slc17a7", "Slc17a6", "Camk2a", "Satb2", "Neurod6", "Tbr1", "Grin1"],
        "Neurons_Inhibitory": ["Gad1", "Gad2", "Slc32a1", "Pvalb", "Sst", "Vip", "Npy"],
        "MSN_D1": ["Drd1", "Tac1", "Pdyn", "Chrm4"],
        "MSN_D2": ["Drd2", "Adora2a", "Penk", "Gpr6"],
        "Astrocytes": ["Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "Slc1a2", "S100b", "Sox9"],
        "Oligodendrocytes": ["Mbp", "Plp1", "Mog", "Olig2", "Olig1", "Cnp", "Mag"],
        "Microglia": ["Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Itgam", "Csf1r", "Hexb"],
        "Endothelial": ["Cldn5", "Pecam1", "Flt1", "Tie1", "Vwf", "Cdh5"],
        "Ependymal": ["Foxj1", "Dnah5", "Dnah9", "Dnah11", "Ccdc153", "Rsph1"],
        "OPC": ["Pdgfra", "Cspg4", "Sox10", "Gpr17"],
    }

    R_TEMPLATE = '''#!/usr/bin/env Rscript
library(DESeq2)
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
    genotype = factor(ifelse(grepl("^{control_prefix}_", sample_names), "{control_label}", "{treatment_label}"),
                      levels=c("{control_label}","{treatment_label}")),
    row.names = sample_names
)

dds <- DESeqDataSetFromMatrix(counts_int, meta, ~ genotype)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=TRUE)
expr <- assay(vsd)
rownames(expr) <- gene_names[match(rownames(expr), rownames(counts_int))]
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "",]
expr <- expr[!duplicated(rownames(expr)),]

# Z-score
expr_z <- t(scale(t(expr)))

# Markers
markers <- {markers_r}

cell_scores <- matrix(NA, nrow=length(markers), ncol=ncol(expr_z))
rownames(cell_scores) <- names(markers)
colnames(cell_scores) <- colnames(expr_z)

for (ct in names(markers)) {{
    present <- markers[[ct]][markers[[ct]] %in% rownames(expr_z)]
    cat(sprintf("  %s: %d/%d markers\\n", ct, length(present), length(markers[[ct]])))
    if (length(present) >= 2) cell_scores[ct,] <- colMeans(expr_z[present,,drop=FALSE], na.rm=TRUE)
    else if (length(present) == 1) cell_scores[ct,] <- expr_z[present,]
}}
valid <- !is.na(cell_scores[,1])
cell_scores <- cell_scores[valid,,drop=FALSE]
write.csv(as.data.frame(cell_scores), file.path(out_dir, "cell_type_scores.csv"))

# Heatmap
anno_c <- data.frame(Genotype=meta$genotype, row.names=rownames(meta))
png(file.path(out_dir, "cell_type_heatmap.png"), width=800, height=500)
pheatmap(cell_scores, annotation_col=anno_c, fontsize_row=10, main="Cell Type Marker Scores")
dev.off()

# Test WT vs KO
ct_tests <- data.frame()
for (ct in rownames(cell_scores)) {{
    tt <- tryCatch(t.test(cell_scores[ct, meta$genotype=="{treatment_label}"],
                          cell_scores[ct, meta$genotype=="{control_label}"]),
                   error=function(e) NULL)
    if (!is.null(tt))
        ct_tests <- rbind(ct_tests, data.frame(CellType=ct,
            mean_ctrl=round(mean(cell_scores[ct, meta$genotype=="{control_label}"]),4),
            mean_trt=round(mean(cell_scores[ct, meta$genotype=="{treatment_label}"]),4),
            diff=round(tt$estimate[1]-tt$estimate[2],4), pvalue=signif(tt$p.value,3)))
}}
ct_tests <- ct_tests[order(ct_tests$pvalue),]
write.csv(ct_tests, file.path(out_dir, "cell_type_tests.csv"), row.names=FALSE)
print(ct_tests)

# Boxplot
ct_melt <- melt(cell_scores)
colnames(ct_melt) <- c("CellType","Sample","Score")
ct_melt$Genotype <- meta$genotype[match(ct_melt$Sample, rownames(meta))]
p <- ggplot(ct_melt, aes(x=CellType, y=Score, fill=Genotype)) +
    geom_boxplot(alpha=0.7) + geom_jitter(width=0.1, size=1.5, alpha=0.7) +
    scale_fill_manual(values=c("{control_label}"="#377EB8","{treatment_label}"="#E41A1C")) +
    labs(title="Cell Type Composition", y="Marker Score") +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_dir, "cell_type_boxplot.png"), p, width=10, height=6, dpi=150)

cat("Deconvolution complete.\\n")
'''

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_script(self, counts_file: str, markers: dict = None,
                        control_prefix: str = "WT", control_label: str = "WT",
                        treatment_label: str = "KO") -> str:
        if markers is None:
            markers = self.BRAIN_MARKERS

        # Convert to R list format
        markers_r = "list(\n" + ",\n".join(
            f'    {k} = c({", ".join(repr(g) for g in v)})'
            for k, v in markers.items()
        ) + "\n)"

        script = self.R_TEMPLATE.format(
            output_dir=str(self.output_dir),
            counts_file=counts_file,
            markers_r=markers_r,
            control_prefix=control_prefix,
            control_label=control_label,
            treatment_label=treatment_label
        )
        path = self.output_dir / "run_deconvolution.R"
        path.write_text(script)
        return str(path)
