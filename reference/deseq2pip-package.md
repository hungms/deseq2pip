# DESeq2 Pipeline for RNA-seq Analysis

A comprehensive R package for RNA-seq data analysis using DESeq2. This
package provides a streamlined workflow for quality control,
differential expression analysis, and gene set enrichment analysis.

## Main Functions

- `run_deseq2_pip()`:

  Main pipeline function that runs the complete analysis workflow

- [`run_diffexp()`](https://hungms.github.io/deseq2pip/reference/run_diffexp.md):

  Performs differential expression analysis between two groups

- [`run_gsea()`](https://hungms.github.io/deseq2pip/reference/run_gsea.md):

  Performs gene set enrichment analysis

- [`plot_volcano()`](https://hungms.github.io/deseq2pip/reference/plot_volcano.md):

  Creates volcano plots for differential expression results

- `plot_gene_expression()`:

  Visualizes gene expression patterns

## Quality Control Functions

- [`remove_xy_genes()`](https://hungms.github.io/deseq2pip/reference/remove_xy_genes.md):

  Removes genes on X and Y chromosomes

- [`remove_mt_genes()`](https://hungms.github.io/deseq2pip/reference/remove_mt_genes.md):

  Removes mitochondrial genes

- [`remove_low_expression()`](https://hungms.github.io/deseq2pip/reference/remove_low_expression.md):

  Filters out lowly expressed genes

- [`check_library()`](https://hungms.github.io/deseq2pip/reference/check_library.md):

  Generates library size distribution plots

- [`run_pca()`](https://hungms.github.io/deseq2pip/reference/run_pca.md):

  Performs principal component analysis

- [`run_distance()`](https://hungms.github.io/deseq2pip/reference/run_distance.md):

  Calculates and visualizes sample distances

## Data Formatting Functions

- `import_nfcore_dds()`:

  Imports DESeq2 object from nfcore pipeline

- `import_msigdb()`:

  Imports MSigDB gene sets

- `summarize_genes()`:

  Aggregates isoform-level expression to gene-level

- `format_enrichmentmap()`:

  Formats results for EnrichmentMap visualization

## Dependencies

The package depends on the following R packages:

- DESeq2

- biomaRt

- clusterProfiler

- ggplot2

- pheatmap

- msigdbr

- fgsea

- ggpubr

## Usage
