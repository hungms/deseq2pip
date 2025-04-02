# deseq2pip 0.1.0

## New Features

### Core Pipeline Functions
- Added `run_diffexp_pip()` for running the complete differential expression analysis pipeline
- Added `run_qc_pip()` for running the quality control pipeline
- Added `run_gsea_pip()` for running the GSEA analysis pipeline
- Added `run_deg_pip()` for running the complete DEG analysis pipeline

### Quality Control Functions
- Added `remove_xy_genes()` for removing X and Y chromosome genes
- Added `remove_mt_genes()` for removing mitochondrial genes
- Added `remove_low_expression()` for filtering lowly expressed genes
- Added `check_library()` for checking library size distribution
- Added `run_pca()` for performing PCA analysis
- Added `run_distance()` for calculating sample distances

### Plotting Functions
- Added `plot_volcano()` for creating volcano plots
- Added `plot_gene_exprs()` for plotting gene expression boxplots
- Added `plot_gsea_barplot()` for visualizing GSEA results

### Utility Functions
- Added `save_tsv()` for saving data frames to TSV files
- Added `save_expression()` for saving expression data
- Added `save_plot()` for saving plots to PDF files
- Added `theme_border()`, `theme_text()`, `theme_gridlines()`, and `facet_aes()` for consistent plot styling

### GSEA Database
- Added support for MSigDB gene sets including:
  - HALLMARK
  - GOBP (GO Biological Process)
  - KEGG
  - REACTOME
  - BIOCARTA
  - TFT (Transcription Factor Targets)
- Support for both human and mouse organisms

## Improvements
- Added customizable color schemes for plots
- Added statistical significance indicators for gene expression plots
- Added support for custom gene set collections in GSEA analysis
- Improved plot aesthetics and formatting
- Added comprehensive documentation and examples

## Bug Fixes
- Fixed issues with file path handling in various functions
- Fixed issues with statistical significance calculations
- Fixed issues with plot dimensions and scaling

## Documentation
- Added comprehensive documentation for all functions
- Added examples for common use cases
- Added detailed parameter descriptions
- Added vignettes demonstrating package usage

## Dependencies
- Requires R >= 4.0.0
- Main dependencies:
  - DESeq2
  - ggplot2
  - dplyr
  - tidyr
  - msigdbr
  - fgsea
  - ggrepel
  - scales
  - rstatix
