FROM rocker/rstudio:latest

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev libxml2-dev libcurl4-openssl-dev libssl-dev \
    libpng-dev libjpeg-dev libfontconfig1-dev libfreetype6-dev \
    libharfbuzz-dev libfribidi-dev libtiff5-dev \
    python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# R packages
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org"); \
    BiocManager::install(c( \
        "DESeq2", "edgeR", "limma", "tximport", "apeglm", "IHW", \
        "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db", "GO.db", \
        "fgsea", "msigdbr", "WGCNA", \
        "EnhancedVolcano", "pheatmap", "RColorBrewer", \
        "ggplot2", "ggrepel", "VennDiagram", "openxlsx", "reshape2" \
    ), ask=FALSE, update=FALSE)'

# Python packages
RUN pip3 install anthropic pandas numpy jinja2 pyyaml rich click

WORKDIR /data
CMD ["bash"]
