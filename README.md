# deseq2pip: A DESeq2 Pipeline for RNA/ATAC-seq Analysis

[![R-CMD-check](https://github.com/hungms/deseq2pip/workflows/R-CMD-check/badge.svg)](https://github.com/hungms/deseq2pip/actions)
[![pkgdown](https://github.com/hungms/deseq2pip/workflows/pkgdown/badge.svg)](https://github.com/hungms/deseq2pip/actions)

## Overview

`deseq2pip` is a comprehensive R package that streamlines RNA-seq data analysis by combining DESeq2-based differential expression analysis with downstream functional analysis and visualization. The package provides a modular yet integrated workflow for quality control, differential expression analysis, gene set enrichment analysis, and visualization of results.

## Documentation

For detailed usage and documentation of all functions, please visit our latest [documentation](https://hungms.github.io/deseq2pip).

## Features

- **Modular Pipeline**: Separate analysis steps that can be run independently or as a complete workflow
  - Quality control and data preparation
  - Differential expression analysis
  - Functional enrichment analysis
  - Visualization of results

- **Data Preprocessing**:
  - Filtering of lowly expressed genes
  - Options to remove sex chromosome genes and mitochondrial genes
  - Quality control plots (PCA, sample distance)

- **Differential Expression Analysis**:
  - Automated pairwise comparisons between experimental groups
  - Integrated DESeq2 workflow with convenient parameter settings
  - Comprehensive result tables with gene-level functional annotation

- **Functional Analysis**:
  - Gene set enrichment analysis (GSEA) using MSigDB gene sets
  - Support for both human and mouse organisms
  - Customizable ranking metrics and significance thresholds

- **Visualization**:
  - Publication-ready volcano plots
  - Customizable gene expression plots
  - GSEA barplots for pathway visualization
  - Support for output formatting for Cytoscape EnrichmentMap


## Workflow
![deseq2pip workflow](deseq2pip_workflow.png)

## Installation

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install deseq2pip from GitHub
devtools::install_github("hungms/deseq2pip")
```

## Dependencies

The package relies on several Bioconductor packages:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages
BiocManager::install(c("DESeq2", "clusterProfiler", "msigdbr", "biomaRt"))
```

## Quick Start

```r
library(deseq2pip)

# load example DESeq2 object
rdata <- system.file("data", "GSE189410.Rdata", package = "deseq2pip")
tx2gene <- system.file("data", "GSE189410_tx2gene.tsv", package = "deseq2pip")
dds <- import_nfcore_rna(rdata = rdata, tx2gene = tx2gene)

# set group variable to compare between, with the first level set as the reference
dds$Group2 <- factor(dds$Group2, c('IgM', 'IgG', 'IgA'))

# set dds design with desired variables
design(dds) <- ~ Group2

# Set directory to store results
save_dir <- "/Users/hungm/Documents/development/deseq2pip/tests/pipeline/"

# Run the complete pipeline
run_deseq2_pip(
    dds,
    assaytype = "RNA",
    experiment = "GSE189410",
    remove_xy = TRUE,
    remove_mt = TRUE,
    org = "mouse",
    group_by = "Group2",
    save_dir = save_dir)

# See output directory structure
list.dir(save_dir)

```

## License

This package is distributed under the MIT License.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

For questions or issues, please open an issue on GitHub.