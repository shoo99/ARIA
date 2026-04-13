"""
PPI Network Module

Protein-protein interaction network analysis using STRING DB API.
"""

from pathlib import Path


class PPIModule:
    """PPI network analysis via STRING DB."""

    R_TEMPLATE = '''#!/usr/bin/env Rscript
library(ggplot2)
library(ggrepel)

out_dir <- "{output_dir}"
de_file <- "{de_results_file}"

de <- read.csv(de_file)
sig <- subset(de, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > {lfc_cutoff})
genes <- na.omit(sig$gene_name)
cat("DEGs for PPI:", length(genes), "\\n")

if (length(genes) < 5) {{ cat("Too few genes for PPI\\n"); quit(save="no") }}

genes_str <- paste(genes, collapse="%0d")
id_url <- paste0("https://string-db.org/api/tsv/get_string_ids?identifiers=",
                  genes_str, "&species={taxid}&limit=1")
id_data <- tryCatch(read.delim(url(id_url), stringsAsFactors=FALSE), error=function(e) NULL)

if (is.null(id_data) || nrow(id_data)==0) {{ cat("STRING ID mapping failed\\n"); quit(save="no") }}

string_ids <- paste(id_data$stringId, collapse="%0d")

# Network interactions
net_url <- paste0("https://string-db.org/api/tsv/network?identifiers=", string_ids,
                  "&species={taxid}&required_score=400")
net <- tryCatch(read.delim(url(net_url), stringsAsFactors=FALSE), error=function(e) NULL)

if (!is.null(net) && nrow(net) > 0) {{
    cat("Interactions:", nrow(net), "\\n")
    write.csv(net[,c("preferredName_A","preferredName_B","score")],
              file.path(out_dir, "PPI_interactions.csv"), row.names=FALSE)

    # Hub genes
    all_nodes <- unique(c(net$preferredName_A, net$preferredName_B))
    degree <- sapply(all_nodes, function(g) sum(net$preferredName_A==g | net$preferredName_B==g))
    hub <- data.frame(gene=names(degree), degree=as.numeric(degree))
    hub$log2FC <- sig$log2FoldChange[match(hub$gene, sig$gene_name)]
    hub <- hub[order(-hub$degree),]
    write.csv(hub, file.path(out_dir, "PPI_hub_genes.csv"), row.names=FALSE)

    hub_plot <- hub[!is.na(hub$log2FC),]
    if (nrow(hub_plot) > 0) {{
        p <- ggplot(hub_plot, aes(x=log2FC, y=degree, label=gene)) +
            geom_point(aes(color=log2FC>0, size=degree), alpha=0.7) +
            geom_text_repel(data=head(hub_plot,15), size=3, max.overlaps=20) +
            scale_color_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("Down","Up")) +
            labs(title="PPI Hub Genes", x="log2FC", y="Degree") + theme_bw()
        ggsave(file.path(out_dir, "PPI_hub_genes.png"), p, width=9, height=6, dpi=150)
    }}

    # Network image
    img_url <- paste0("https://string-db.org/api/image/network?identifiers=", string_ids,
                      "&species={taxid}&required_score=400&network_flavor=confidence")
    tryCatch(download.file(img_url, file.path(out_dir, "PPI_network.png"), mode="wb", quiet=TRUE),
             error=function(e) NULL)
}}

cat("PPI analysis complete.\\n")
'''

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_script(self, de_results_file: str, taxid: int = 10090,
                        lfc_cutoff: float = 0.5) -> str:
        script = self.R_TEMPLATE.format(
            output_dir=str(self.output_dir),
            de_results_file=de_results_file,
            taxid=taxid,
            lfc_cutoff=lfc_cutoff
        )
        path = self.output_dir / "run_PPI.R"
        path.write_text(script)
        return str(path)
