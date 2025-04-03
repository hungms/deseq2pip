#' DESeq2 Pipeline for RNA-seq Analysis
#'
#' A comprehensive R package for RNA-seq data analysis using DESeq2. This package
#' provides a streamlined workflow for quality control, differential expression analysis,
#' and gene set enrichment analysis.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{run_deseq2_pip()}}{Main pipeline function that runs the complete analysis workflow}
#'   \item{\code{run_diffexp()}}{Performs differential expression analysis between two groups}
#'   \item{\code{run_gsea()}}{Performs gene set enrichment analysis}
#'   \item{\code{plot_volcano()}}{Creates volcano plots for differential expression results}
#'   \item{\code{plot_gene_expression()}}{Visualizes gene expression patterns}
#' }
#'
#' @section Quality Control Functions:
#' \describe{
#'   \item{\code{remove_xy_genes()}}{Removes genes on X and Y chromosomes}
#'   \item{\code{remove_mt_genes()}}{Removes mitochondrial genes}
#'   \item{\code{remove_low_expression()}}{Filters out lowly expressed genes}
#'   \item{\code{check_library()}}{Generates library size distribution plots}
#'   \item{\code{run_pca()}}{Performs principal component analysis}
#'   \item{\code{run_distance()}}{Calculates and visualizes sample distances}
#' }
#'
#' @section Data Formatting Functions:
#' \describe{
#'   \item{\code{import_nfcore_dds()}}{Imports DESeq2 object from nfcore pipeline}
#'   \item{\code{import_msigdb()}}{Imports MSigDB gene sets}
#'   \item{\code{summarize_genes()}}{Aggregates isoform-level expression to gene-level}
#'   \item{\code{format_enrichmentmap()}}{Formats results for EnrichmentMap visualization}
#' }
#'
#' @section Dependencies:
#' The package depends on the following R packages:
#' \itemize{
#'   \item DESeq2
#'   \item biomaRt
#'   \item clusterProfiler
#'   \item ggplot2
#'   \item pheatmap
#'   \item msigdbr
#'   \item fgsea
#'   \item ggpubr
#' }
#'
#' @section Usage:
#' \dontrun{
#' # Load the package
#' library(deseq2pip)
#'
#' # Run the complete pipeline
#' dds <- run_deseq2_pip(dds, experiment = "my_experiment")
#'
#' # Run individual analyses
#' res <- run_diffexp(dds)
#' gsea_results <- run_gsea(res)
#' }
#'
#' @name deseq2pip-package
NULL
