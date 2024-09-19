---
title: "Data Wrangling CHD Cohort README"
author: "Jacqueline Penaloza and Peter White"
date: "`r Sys.Date()`"
department: "The Office of Data Sciences"
organization: "Abigail Wexner Research Institute, Nationwide Children's Hospital"
address: "575 Children's Crossroad, Columbus, OH 43215"
output:
  rmdformats::downcute:
    code_folding: show           # Show code by default with folding
    self_contained: true         # Ensure the HTML file is self-contained
    thumbnails: false            # Disable thumbnails
    lightbox: true               # Enable lightbox for images
    gallery: false               # Disable gallery view for images (replaces thumbnails)
    downcute_theme: "default"    # Use the light theme (options: “default”, “chaos”)
    number_sections: true        # Enable numbered sections
    toc_depth: 3                 # Set the depth of the table of contents
    css: "Penaloza_Custom_Styling.css"  # Link to the external CSS file
editor_options:
  markdown:
    wrap: sentence               # Wrap lines at sentence boundaries
  chunk_output_type: console     # Show code chunk output in the console
pkgdown:
  as_is: true                    # Maintain formatting with pkgdown
---

# CCVM Long Non-Coding RNA (lncRNA) Analysis

This repository contains R scripts for the analysis of long non-coding RNAs (lncRNAs) associated with copy number variants (CNVs) in congenital heart disease (CHD) patients from the CCVM cohort. The analysis includes transcript quantification, network construction using WGCNA, and functional enrichment analysis of lncRNAs involved in CHD.

## Project Overview

In this project, we aim to investigate the role of lncRNAs impacted by CNVs in CHD. Using the WGCNA method, we identified key gene modules and candidate lncRNAs associated with heart development. This repository contains R Markdown notebooks that guide the data processing, analysis, and visualization steps:

1. **Penaloza_CCVM_Notebook_1A.Rmd**: Preprocessing and normalization of gene expression data (TPM), including filtering by expression threshold and initial clustering of genes for network analysis.
   
2. **Penaloza_CCVM_Notebook_1B.Rmd**: Construction of weighted gene co-expression networks (WGCNA), identification of heart-specific modules, and functional enrichment analysis of candidate lncRNAs.

## Repository Contents

- **Penaloza_CCVM_Notebook_1A.Rmd**: Contains code for the preprocessing of RNA-Seq data, thresholding TPM values, and data normalization.
- **Penaloza_CCVM_Notebook_1B.Rmd**: Contains code for WGCNA analysis, identification of hub genes, and visualization of key gene modules.
- **InputData**: This folder contains the CCVM data file from which all analysis is dependent. `02032023_CCVM_abnormalEcho.Rds` is a binary R file required for the script to run. `02032023_CCVM_abnormalEcho.xlsx` is the data in Microsoft Excel format.
- **/PublicData**: This folder contains all of the publically available data that is required for the analysis. This includes:
-- *lncExpDB_E-MTABGeneTPM.tsv*: Normalized TPM values for the developmental time series data
-- *E-MTAB-6814.sdrf.txt*: Metadata for the developmental time series
-- *LncBookv1.9_GENCODEv33_GRCh38.gtf.gz*: GTF file used for RNA-Seq transcript quantification
-- *UCSC_hg19ToHg38.over.chain*: Chain file to liftover genomic coordinates from GRCh37 to GRCh38
-- *LncBook_id_conversion.csv*: LncBook provided conversions of LncBook accessions
-- *RNACentral_lncbook.tsv*: RNACentral to LncBook accessions
-- *RNACentral_ensembl.tsv*: RNACenral to Ensembl accessions
-- *chdgene_table.csv*: List of 142 known CHD genes from CHDgene website
-- *Cotney_CircRes_316709_online_table_v.xlsx*: VanOudenhove et al., 2020. Table V
- **/Results**: This folder will contain output files such as gene module lists, hub gene lists, and GO enrichment analysis results after running the scripts.
- **/Figures**: This folder will contain plots and visualizations generated during the analysis. All are found in our manuscript.
- **/


## Dependencies

The following R packages are required to run the analysis:

### Data Manipulation and Tidy Data

- `tidyr`
- `dplyr`
- `stringr`
- `purrr`
- `tibble`
- `readr`
- `data.table`
- `reshape2`
- `foreach`
- `doParallel`

### Visualization and Plotting

- `ggplot2`
- `grid`
- `gridExtra`
- `ggdendro`
- `ggrepel`
- `Cairo`
- `viridis`
- `gplots`
- `ggpubr`
- `circlize`
- `eulerr`

### Genomic Data and Network Analysis

- `GenomicRanges`
- `biomaRt`
- `rtracklayer`
- `WGCNA`
- `RCy3`
- `igraph`
- `clusterProfiler`
- `GOSemSim`
- `org.Hs.eg.db`
- `gprofiler2`

### File Handling and Export

- `readxl`
- `writexl`
- `R.utils`
- `htmlwidgets`
- `webshot2`
- `chromote`

### Reporting and Formatting

- `knitr`
- `kableExtra`
- `rmdformats`
- `skimr`
- `here`

Packages available on CRAN are installed using `install.packages()`, and those available on Bioconductor are installed using `BiocManager::install()`.

### CRAN Packages:

```r
# Install CRAN packages
install.packages(c(
  "tidyr", "dplyr", "stringr", "purrr", "tibble", "readr", "data.table",
  "reshape2", "foreach", "doParallel", "ggplot2", "grid", "gridExtra",
  "ggdendro", "ggrepel", "Cairo", "viridis", "gplots", "ggpubr", "circlize",
  "eulerr", "readxl", "writexl", "R.utils", "htmlwidgets", "webshot2",
  "chromote", "knitr", "kableExtra", "rmdformats", "skimr", "here"
))
```

### Bioconductor Packages:

```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges", "biomaRt", "rtracklayer", "WGCNA", "RCy3", "clusterProfiler",
  "GOSemSim", "org.Hs.eg.db", "gprofiler2"
))
```

### Explanation:

- **CRAN Packages**: These packages are installed using `install.packages()`. Simply run the first command in R to install all the CRAN packages in one go.
- **Bioconductor Packages**: Bioconductor packages are installed using `BiocManager::install()`. Make sure the `BiocManager` package is installed first, and then install the Bioconductor-specific packages listed.

This should install all the necessary dependencies for the CNV-lncRNA analysis project.

### Setup

Ensure that the data files required for analysis (e.g., TPM matrices, metadata) are available in the appropriate directories specified in the scripts.

## Usage

1. **Penaloza_CCVM_Notebook_1A.Rmd**: Preprocessing and Data Normalization
   - Run this notebook to preprocess RNA-Seq data, normalize the TPM values, and threshold the data for WGCNA analysis.
   
   Command (in RStudio or terminal):
   ```r
   rmarkdown::render("Penaloza_CCVM_Notebook_1A.Rmd")
   ```

2. **Penaloza_CCVM_Notebook_1B.Rmd**: WGCNA and Hub Gene Identification
   - Run this notebook to perform WGCNA, identify heart-specific gene modules, and find hub genes.

   Command (in RStudio or terminal):
   ```r
   rmarkdown::render("Penaloza_CCVM_Notebook_1B.Rmd")
   ```

## Output

- **Gene Modules**: Lists of genes associated with each WGCNA module, saved in `/Results`.
- **Hub Genes**: Identified hub genes within key heart-related modules.
- **GO Enrichment Analysis**: Functional enrichment results for modules, including GO biological processes, molecular functions, and cellular components.
- **Figures**: Key plots such as gene module dendrograms, eigengene expression patterns, and co-expression networks (saved in `/Figures`).
- **Results**: The results of the filtering and liftover process from Notebook 1A, `Penaloza_lncRNA_CCVM_hg38_liftover.rds` and the intersection of this data with expressed lncRNAs, `Penaloza_expressed_lncRNA_CCVM_hg38_liftover.rds`. The results of processing the LncBook GTF (`NCH_LncBookv1.9_GENCODEv33_GRCh38.gtf.gz` and `NCH_LncBookv1.9_GENCODEv33_GRCh38.gff.gz`). Supplemental Data Tables and results tables (saved in `/Results`).

## Network Visualization

The function `create_and_view_network()` within the notebook `Penaloza_CCVM_Notebook_1B.Rmd` allows for the creation and visualization of gene co-expression networks. The network can be viewed using Cytoscape, and the function saves the network as a GraphML file for further analysis.
# HeartiLNC
