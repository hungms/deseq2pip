# NEWS
## v2.1.1 (2026-04-24)

### Bug Fixes

- **`run_diffexp()` one-vs-rest (pooled) comparisons:**
  - Merging factor levels for comparisons containing `+` now uses exact level replacement instead of regular-expression `gsub()`, avoiding mis-mapping when level names are substrings of each other.
  - Before calling `DESeq()`, the pipeline assigns the design to the temporary `DESeqDataSet` and checks rank using the same `stats::model.matrix(design(object), ...)` construction as DESeq2’s `designAndArgChecker` (with `betaPrior = FALSE`). Pooled comparisons also error if any sample’s `var` level is not listed on either side of the contrast after merging.
  - **Shared `dds` mutation:** each call to `run_diffexp()` now replaces only `colData(temp_dds)[[var]]` with a duplicated factor so pooled one-vs-rest cannot collapse levels on the parent `dds` reused by `run_diffexp_pip()` (which broke later pairwise contrasts / `cleanContrast`). A full `colData()` replacement was avoided because it strips `mcols(colData)` and breaks DESeq2 outlier handling.

## v2.1.0 (2026-03-17)
### New Features
- **`comparisons` argument in main pipeline functions:**
  - `run_rna_pip()` and `run_atac_pip()` now accept an explicit `comparisons` argument, allowing users to specify a custom subset of comparisons to run instead of the default auto-generated full set.
  - When `comparisons = NULL` (default), all pairwise and one-to-all comparisons are automatically generated from the factor levels of `var`.
  - Supports standard pairwise comparisons (e.g. `"B_vs_A"`), combined-group comparisons (e.g. `"IgA+IgG_vs_IgM"`), or any mix of both.

### Bug Fixes
- **Fixed `validate_comparison()` for combined-group comparisons:**
  - `validate_comparison()` now correctly handles comparisons containing `+` (e.g. `"IgA+IgG_vs_IgM"`) by splitting each side of `_vs_` further on `+` before checking group membership. Previously, combined-group comparisons that passed `validate_comparisons()` would fail at the inner `validate_comparison()` call inside `run_diffexp()`.

## v1.0.4 (2025-12-23)
### Minor Changes
- **Updated `generate_comparisons` function:**
  - set comparisons with regards to factor levels

## v1.0.3 (2025-01-15)
### Major Changes
- **Updated `remove_low_expression()` function:**
  - Changed parameter from `quantile` to `min_count` for more intuitive filtering
  - Now filters based on raw count values (default: min_count = 10) instead of quantile-based filtering on transformed data
  - Added percentage reporting of genes remaining after filtering
  - Enhanced visualization with separate before/after filtering density plots
  - Improved plot naming: "low_expression_before.pdf" and "low_expression_after.pdf"

- **Updated validation functions:**
  - Added `validate_min_count()` function for validating minimum count thresholds
  - Removed `validate_quantile()` function (no longer needed)

- **Updated pipeline functions:**
  - `run_qc_pip()`, `run_rna_pip()`, and `run_atac_pip()` now use `min_count` parameter
  - All functions consistently use the new parameter naming convention

- **Updated documentation:**
  - Updated vignettes and examples to reflect new `min_count` parameter
  - Improved function documentation and parameter descriptions

## v1.0.2 (2025-06-25)
Change all plotr dependencies to ggprism.

## v1.0.1 (2025-05-30)
reverse groups in all one-to-all comparisons.

## v1.0.0 (2025-05-20)
### Bugfix
- Make compatible with strpip v0.1.3


## v1.0.0 (2025-04-30)

### Major Refactor and Pipeline Redesign
- **Unified and modularized pipeline:**
  - `run_rna_pip()` and `run_atac_pip()` have been refactored for clarity, modularity, and robust argument validation. All input checks are now handled by dedicated `validate_*` subfunctions.
  - Arguments `org` and `group_by` are now required and must be explicitly specified.
  - The `one_to_all` argument is now handled internally; one-to-all and pairwise comparisons are automatically determined based on group levels.
  - **Logging is now always enabled by default:** all output, messages, warnings, and errors are captured for every run of the main pipeline functions, with no need for a user toggle.
  - Removed legacy files: `R/pipeline_functions.R`, `R/format_functions.R`, and `R/quality_control.R`.

### Quality Control Functions

- **`run_qc_pip()`**
  - Unified and streamlined quality control pipeline for both RNA-seq and ATAC-seq workflows.
  - Handles gene/peak filtering, quality control plots, and expression data saving in a single function.
  - Now supports explicit argument validation for organism, group, and quantile threshold.
  - Automatically removes X/Y and mitochondrial genes if requested.
  - Improved low-expression filtering with customizable quantile threshold.
  - Saves all relevant QC outputs (filtered data, plots, and summary statistics) in organized directories.

- **`remove_xy_genes()` / `remove_mt_genes()`**
  - Enhanced logic for robust removal of X, Y, and mitochondrial genes based on organism.
  - Improved error handling for missing or malformed gene annotations.

- **`remove_low_expression()`**
  - Improved quantile-based filtering for lowly expressed genes.
  - Now provides informative messages and summary statistics on filtering results.

- **`check_library()`**
  - Enhanced library size distribution checks and reporting.
  - Improved error messages for outlier detection.

- **`save_expression()`**
  - Now saves DESeq2 object, raw counts, normalized expression, and class labels in a consistent and organized manner.
  - Improved file naming and directory structure for easier downstream analysis.

- **`run_dist_pip()`**
  - Unified distance and PCA analysis for both RNA-seq and ATAC-seq.
  - Supports batch correction and custom color palettes.
  - Improved plot aesthetics and output organization.

- **Validation and Error Handling**
  - All QC functions now use centralized `validate_*` helpers for argument and data validation.
  - Improved error messages and user feedback throughout the QC pipeline.

- **Logging and Reproducibility**
  - All QC steps are now logged by default, including function calls, arguments, and all messages/warnings/errors, for full reproducibility.

### Differential Expression Functions
- **`run_diffexp_pip()`**
  - Now supports both pairwise and one-to-all comparisons in a unified interface.
  - Uses new `generate_comparisons()` to robustly generate all group comparisons.
  - Differential expression results are saved in structured subdirectories for each comparison.
  - MA and volcano plots are generated and saved for each comparison automatically.
  - Wrapper function `run_diffexp_wrapper()` added for single-comparison analysis and plotting.
  - Improved merging of row metadata and annotation.
  - Enhanced error handling and validation for group and comparison arguments.
- **`run_diffexp()`**
  - Now requires explicit `org`, `group_by`, and `comparison` arguments.
  - Subsets DESeq2 object to the correct samples for each comparison.
  - Improved annotation and filtering of results.
- **`read_diffexp()`**
  - Replaces `read_diffexp_list()`. Reads and merges results from all comparisons in a directory.

### GSEA Functions
- **`run_gsea_pip()` / `run_gsea_wrapper()` / `run_gsea()`**
  - GSEA pipeline now modularized: runs GSEA for all comparisons and gene set collections.
  - New `import_msigdbr()` for importing MSigDB gene sets for human or mouse.
  - GSEA results are saved in structured subdirectories, with both RDS and TSV outputs.
  - Barplots for GSEA results are generated for each collection and comparison.
  - Improved validation for gene set input and result structure.
- **`read_gsea_rds()` / `read_gsea_tsv()`**
  - Replaces `read_gsea_rds_list()` and `read_gsea_tsv_list()`. Reads and merges GSEA results from all comparisons and collections.

### Plotting Functions
- **MA and Volcano Plots**
  - `plot_ma()` and `plot_volcano()` are now internal to the pipeline and automatically called for each comparison.
  - Improved aesthetics, labeling, and filtering for large datasets.
- **GSEA Barplots**
  - `plot_gsea_barplot()` now supports filtering for significant gene sets and improved labeling.
- **Gene Expression Plots**
  - `plot_gene_exprs()` now requires explicit `group_by` and supports custom color palettes via `pal` argument.
- **ATAC-seq Annotation Plots**
  - `plot_atac_annot()` and related list functions are now internal and called automatically in ATAC-seq pipelines.

### Logging and Reproducibility
- **`log_renv()`**
  - Captures and logs the renv environment snapshot for each run.
- **`log_output()`**
  - Captures and logs all output, messages, warnings, and errors from pipeline runs.
  - Logging is always enabled for all major pipeline functions.

### Enhancements and Validation
- Centralized all argument validation in new `validate_*` helper functions for consistency and robustness.
- Improved error messages and user feedback throughout the pipeline.
- Improved directory creation and file path handling for all save and log functions.
- Updated vignettes and documentation to reflect new argument names, usage, and workflow patterns.

### Breaking Changes
- **Function signatures:**
  - Main pipeline functions (`run_rna_pip`, `run_atac_pip`, `run_diffexp`, etc.) now require explicit `org` and `group_by` arguments.
  - Some arguments have been renamed, reordered, or removed for clarity and consistency.
- **Removed deprecated/legacy functions and files:**
  - `R/pipeline_functions.R`, `R/format_functions.R`, `R/quality_control.R` and their man pages have been deleted.
  - Many old plotting and utility functions are now internal and not exported.

### Bug Fixes
- Fixed issues with directory creation and file path handling in all save and logging functions.
- Improved robustness of GSEA and annotation plotting for edge cases.
- Fixed merging of row metadata in differential expression results.
- Improved error handling for missing or invalid arguments in all main functions.

### Documentation
- Updated vignettes and function documentation for all major changes.
- Improved examples for both RNA-seq and ATAC-seq workflows.
- All function documentation now reflects new argument names and usage patterns.

## v0.1.3
### New Features
- Split pipeline functionality into specialized `run_rna_pip()` and `run_atac_pip()` functions for better workflow separation and clarity
- Added `plot_atac_annot()` function for visualizing genomic annotations of ATAC-seq peaks as pie charts, with separate views for all DE peaks, upregulated peaks, and downregulated peaks
- Added support for displaying package version when running the complete pipeline with `packageVersion("deseq2pip")` in run_deseq2_pip()

### Enhancements
- Enable `pals` argument in `run_dist_pip()` function and related functions
- Change `nobatch` to `default` labelling for batch-correction
- Remove all unnecessary instance of `experiment` argument
- Create new `enrichmentmap` output directory for `format_enrichmentmap()` function for convenient data loading onto Cytoscape
- Change default `order` to `pxfc`, set `rank` as a column containing values for either `log2FoldChange`, `padj`, or `pxfc`
- Improved filtering of differentially expressed peaks based on adjustable thresholds (padj < 0.05 and |log2FoldChange| > 0.5 by default)
- Add x=0 & y=0 lines for PCA plot

### Bugfix
- Enable dataframe input for `custom_msigdb` argument in `run_gsea()`

## v0.1.2
### Enhancements
- Change `diffexp_deseq2_wald.tsv`/`diffexp_deseq2_wald_rank.tsv` to `diffexp_DESeq2.tsv`/`diffexp_DESeq2_rank.tsv` for better readability & compatibility.
- Improve `plot_gsea_enriched` title & subtitle specification
- Improve `plot_gsea_barplot` title & subtitle specification

## v0.1.1

### New Features
- Added compatibility with ATAC-seq data analysis with TSS-specific peak detection
- Added optional batch correction methods using limma to correct normalized counts
- Added `run_dist_pip()` function for sample distance analysis and visualization
- Added `plot_ma()` and `plot_gsea_enriched()` to improve gene/pathway visualization
- Added custom color specification (`cols` parameter) for PCA and gene expression plots
- Added automatic LFC shrinkage with `lfcshrink` for more accurate fold change estimates

### Enhancements

#### Core Functionality
- Enhanced `run_diffexp()` to automatically merge with rowData() for better metadata integration
- Improved error handling and fallback options for VST transformation
- Added more informative messages throughout pipeline functions
- Added more robust parameter validation in all pipeline functions
- Added experiment name parameter to most functions for better file organization
- Changed save directory name default from "qualitycontrol" to "qc_results" for consistency

#### ATAC-seq
- Improved handling of ATAC-seq specific workflows and data structures
- Enhanced TSS-peak subsetting for more accurate transcription start site analysis
- Better handling of duplicated gene names in peak analysis

#### Visualization
- Improved MA plots for better visualization
- Improved gene set enrichment plots for better visualization
- Enhanced PCA plots with support for custom colors and improved aesthetics
- Improved gene expression plots with custom color options

### Bug Fixes
- Fixed batch correction implementation to properly handle design formulas
- Fixed incorrect parameter handling in `save_expression()` function
- Fixed duplicated gene name handling in TSS-specific analysis
- Fixed missing parameter validation in pipeline functions
- Fixed Cairo PDF device fallback mechanism
- Fixed a bug in the `remove_xy_gene` and `remove_mt_gene` function that was causing incorrect behavior.

### Documentation
- Added examples for ATAC-seq data analysis
- Added examples for batch correction usage
- Added session information

## v0.1.0

### New Features

#### Core Pipeline Functions
- Added `run_diffexp_pip()` for running the complete differential expression analysis pipeline
- Added `run_qc_pip()` for running the quality control pipeline
- Added `run_deseq2_pip()` for running the complete DEG analysis pipeline
- Added `run_gsea_pip()` for running the GSEA analysis pipeline

#### Quality Control Functions
- Added `remove_xy_genes()` for removing X and Y chromosome genes
- Added `remove_mt_genes()` for removing mitochondrial genes
- Added `remove_low_expression()` for filtering lowly expressed genes
- Added `check_library()` for checking library size distribution
- Added `run_pca()` for performing PCA analysis
- Added `run_distance()` for calculating sample distances

#### Plotting Functions
- Added `plot_volcano()` for creating volcano plots
- Added `plot_gene_exprs()` for plotting gene expression boxplots
- Added `plot_gsea_barplot()` for visualizing GSEA results

#### Utility Functions
- Added `save_tsv()` for saving data frames to TSV files
- Added `save_expression()` for saving expression data
- Added `save_plot()` for saving plots to PDF files
- Added `theme_border()`, `theme_text()`, `theme_gridlines()`, and `facet_aes()` for consistent plot styling

#### GSEA Database
- Added support for MSigDB gene sets including:
  - HALLMARK
  - GOBP (GO Biological Process)
  - KEGG
  - REACTOME
  - BIOCARTA
  - TFT (Transcription Factor Targets)
- Support for both human and mouse organisms

### Enhancements
- Added customizable color schemes for plots
- Added statistical significance indicators for gene expression plots
- Added support for custom gene set collections in GSEA analysis
- Improved plot aesthetics and formatting
- Added comprehensive documentation and examples

### Bug Fixes
- Fixed issues with file path handling in various functions
- Fixed issues with statistical significance calculations
- Fixed issues with plot dimensions and scaling

### Documentation
- Added comprehensive documentation for all functions
- Added examples for common use cases
- Added detailed parameter descriptions
- Added vignettes demonstrating package usage

### Dependencies
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
