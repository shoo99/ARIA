"""
GSEA Module

Gene Set Enrichment Analysis using fgsea with MSigDB gene sets.
Supports Hallmark, GO:BP/CC/MF, KEGG, Reactome.
"""

from pathlib import Path


class GSEAModule:
    """GSEA analysis using fgsea."""

    R_TEMPLATE = '''#!/usr/bin/env Rscript
library(fgsea)
library(msigdbr)
library(ggplot2)

# Load DE results
de_res <- read.csv("{de_results_file}")

# Build ranked list
de_res <- de_res[!is.na(de_res$gene_name) & !is.na(de_res$stat),]
de_res <- de_res[order(-abs(de_res$stat)),]
de_res <- de_res[!duplicated(de_res$gene_name),]
ranked <- sort(setNames(de_res$stat, de_res$gene_name), decreasing=TRUE)
cat("Ranked genes:", length(ranked), "\\n")

# Gene set collections
collections <- list(
    Hallmark = msigdbr(species="{species}", category="H"),
    GO_BP = msigdbr(species="{species}", category="C5", subcategory="GO:BP"),
    GO_CC = msigdbr(species="{species}", category="C5", subcategory="GO:CC"),
    GO_MF = msigdbr(species="{species}", category="C5", subcategory="GO:MF"),
    KEGG = msigdbr(species="{species}", category="C2", subcategory="CP:KEGG_MEDICUS"),
    Reactome = msigdbr(species="{species}", category="C2", subcategory="CP:REACTOME")
)

for (gs_name in names(collections)) {{
    gs_list <- split(collections[[gs_name]]$gene_symbol, collections[[gs_name]]$gs_name)
    set.seed(42)
    res <- fgsea(pathways=gs_list, stats=ranked, minSize=15, maxSize=500, nPermSimple=10000)
    res <- res[order(res$pval),]
    sig_n <- sum(res$padj < 0.05, na.rm=TRUE)
    cat(sprintf("  %s: %d significant (padj<0.05)\\n", gs_name, sig_n))

    # Save
    out <- as.data.frame(res)
    out$leadingEdge <- sapply(out$leadingEdge, paste, collapse=";")
    write.csv(out, sprintf("{output_dir}/GSEA_%s.csv", gs_name), row.names=FALSE)

    # Plot
    if (nrow(res) >= 5) {{
        top <- head(res[order(res$pval),], 25)
        top$pw_clean <- gsub("^HALLMARK_|^GOBP_|^GOCC_|^GOMF_|^KEGG_MEDICUS_|^REACTOME_", "", top$pathway)
        top$pw_clean <- gsub("_", " ", top$pw_clean)
        top$pw_clean <- substr(top$pw_clean, 1, 65)
        p <- ggplot(top, aes(x=reorder(pw_clean, NES), y=NES, fill=padj<0.05)) +
            geom_col() + coord_flip() +
            scale_fill_manual(values=c("TRUE"="#E41A1C", "FALSE"="#999999"), name="padj<0.05") +
            labs(title=paste0("GSEA ", gs_name), x="", y="NES") +
            theme_bw() + theme(axis.text.y=element_text(size=7))
        ggsave(sprintf("{output_dir}/GSEA_%s.png", gs_name), p, width=12, height=8, dpi=150)
    }}
}}

cat("GSEA complete.\\n")
'''

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_script(self, de_results_file: str,
                        species: str = "Mus musculus") -> str:
        """Generate R script for GSEA."""
        script = self.R_TEMPLATE.format(
            de_results_file=de_results_file,
            output_dir=str(self.output_dir),
            species=species
        )
        script_path = self.output_dir / "run_GSEA.R"
        script_path.write_text(script)
        return str(script_path)
