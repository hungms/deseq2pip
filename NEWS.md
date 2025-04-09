# deseq2pip 0.1.2
## Enhancements
- Change `diffexp_deseq2_wald.tsv`/`diffexp_deseq2_wald_rank.tsv` to `diffexp_DESeq2.tsv`/`diffexp_DESeq2_rank.tsv` for better readability & compatibility.
- Improve `run_gsea_enriched()` visualization

# deseq2pip 0.1.1

## New Features
- Added compatibility with ATAC-seq data analysis with TSS-specific peak detection
- Added optional batch correction methods using limma to correct normalized counts
- Added `run_dist_pip()` function for sample distance analysis and visualization
- Added `plot_ma()` and `plot_gsea_enriched()` to improve gene/pathway visualization
- Added custom color specification (`cols` parameter) for PCA and gene expression plots
- Added automatic LFC shrinkage with `lfcshrink` for more accurate fold change estimates

## Enhancements

### Core Functionality
- Enhanced `run_diffexp()` to automatically merge with rowData() for better metadata integration
- Improved error handling and fallback options for VST transformation
- Added more informative messages throughout pipeline functions
- Added more robust parameter validation in all pipeline functions
- Added experiment name parameter to most functions for better file organization
- Changed save directory name default from "qualitycontrol" to "qc_results" for consistency

### ATAC-seq
- Improved handling of ATAC-seq specific workflows and data structures
- Enhanced TSS-peak subsetting for more accurate transcription start site analysis
- Better handling of duplicated gene names in peak analysis

### Visualization
- Improved MA plots for better visualization
- Improved gene set enrichment plots for better visualization
- Enhanced PCA plots with support for custom colors and improved aesthetics
- Improved gene expression plots with custom color options

## Bug Fixes
- Fixed batch correction implementation to properly handle design formulas
- Fixed incorrect parameter handling in `save_expression()` function
- Fixed duplicated gene name handling in TSS-specific analysis
- Fixed missing parameter validation in pipeline functions
- Fixed Cairo PDF device fallback mechanism
- Fixed a bug in the `remove_xy_gene` and `remove_mt_gene` function that was causing incorrect behavior.

## Documentation
- Added examples for ATAC-seq data analysis
- Added examples for batch correction usage
- Added session information

# deseq2pip 0.1.0

## New Features

### Core Pipeline Functions
- Added `run_diffexp_pip()` for running the complete differential expression analysis pipeline
- Added `run_qc_pip()` for running the quality control pipeline
- Added `run_deseq2_pip()` for running the complete DEG analysis pipeline
- Added `run_gsea_pip()` for running the GSEA analysis pipeline

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
