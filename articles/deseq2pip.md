# Get Started

``` r
# load deseq2pip  
suppressPackageStartupMessages(library(deseq2pip))
```

    ## Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
    ## 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

## Run deseq2pip on RNA-seq data

Create a directory for results, currently set to tempdir()

``` r
save_dir <- tempdir()
dir.create(save_dir, recursive = TRUE)
```

    ## Warning in dir.create(save_dir, recursive = TRUE): '/tmp/RtmpuxUX6f' already
    ## exists

Load DESeq2 object processed by nf-core/rnaseq.

``` r
rdata <- system.file("data", "GSE189410.dds.RData", package = "deseq2pip")
tx2gene <- gzfile(system.file("data", "GSE189410.tx2gene.tsv.gz", package = "deseq2pip"))
dds <- import_nfcore_rna(rdata = rdata, tx2gene = tx2gene)
```

Setup dds comparisons and design.

``` r
# set correct group factors
dds$Group2 <- factor(dds$Group2, c('IgM', 'IgG', 'IgA'))
dds$Group3 <- factor(dds$Group3, c(1:4))

# set design
design(dds) <- ~ Group3 + Group2
```

Run RNA pipeline. When `comparisons = NULL` (default), all pairwise and
one-to-all comparisons are automatically generated from the factor
levels of `var`.

``` r
dds <- run_rna_pip(
    dds,
    org = "mouse",
    var = "Group2",
    design = "~ Group3 + Group2",
    batch = "Group3",
    comparisons = NULL,
    remove_xy = TRUE,
    remove_mt = TRUE,
    min_count = 10,
    pals = NULL,
    order = "pxfc",
    save_dir = save_dir)
```

    ## renv not initialized, initializing renv

    ## - The project is out-of-sync -- use `renv::status()` for details.
    ## The following required packages are not installed:
    ## - rmarkdown
    ## Packages must first be installed before renv can snapshot them.
    ## Use `renv::dependencies()` to see where this package is used in your project.
    ## 
    ## The following package(s) will be updated in the lockfile:
    ## 
    ## # CRAN -----------------------------------------------------------------------
    ## - renv   [* -> 1.2.0]
    ## 
    ## The version of R recorded in the lockfile will be updated:
    ## - R      [* -> 4.5.3]
    ## 
    ## - Lockfile written to "/tmp/RtmpuxUX6f/logs/renv.lock".

    ## Running RNA-seq pipeline with DESeq2pip v2.1.0

    ## 
    ## ########################################################
    ## Running Quality Control Pipeline
    ## ########################################################

    ## removing 4050 XY genes out of 54513 total genes...

    ## removing 37 MT genes out of 50463 total genes...

    ## WARNING: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the deseq2pip package.
    ##   Please report the issue at <https://github.com/hungms/deseq2pip/issues>.

    ## WARNING: Removed 10363 rows containing non-finite outside the scale range
    ## (`stat_density()`).
    ## WARNING: Removed 10363 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](deseq2pip_files/figure-html/unnamed-chunk-5-1.png)

    ## after removing low expression genes, 14178 out of 50426 genes remain (28.12%)

![](deseq2pip_files/figure-html/unnamed-chunk-5-2.png)

    ## generating boxplots to check library size distribution...

![](deseq2pip_files/figure-html/unnamed-chunk-5-3.png)

    ## saving DESeq2 object & expressions...

    ## 
    ## ########################################################
    ## Running Distance Analysis Pipeline
    ## ########################################################

    ## performing variance stabilizing transformation (VST)...

    ## performing limma batch correction...

    ## design matrix of interest not specified. Assuming a one-group experiment.

    ## running PCA analysis...

    ## using ntop=500 top features by variance

    ##  - generating PCA plot for Group2...

    ## WARNING: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the deseq2pip package.
    ##   Please report the issue at <https://github.com/hungms/deseq2pip/issues>.

![](deseq2pip_files/figure-html/unnamed-chunk-5-4.png)

    ##  - generating PCA plot for Group3...

![](deseq2pip_files/figure-html/unnamed-chunk-5-5.png)

    ##  - generating PCA plot for Group1...

![](deseq2pip_files/figure-html/unnamed-chunk-5-6.png)

    ## saving PCA results to TSV...

    ## calculating euclidean distance between samples...

    ## generating sample distance heatmap...

    ## saving distance matrix to TSV...

![](deseq2pip_files/figure-html/unnamed-chunk-5-7.png)

    ## running PCA analysis...

    ## using ntop=500 top features by variance

    ##  - generating PCA plot for Group2...

![](deseq2pip_files/figure-html/unnamed-chunk-5-8.png)

    ##  - generating PCA plot for Group3...

![](deseq2pip_files/figure-html/unnamed-chunk-5-9.png)

    ##  - generating PCA plot for Group1...

![](deseq2pip_files/figure-html/unnamed-chunk-5-10.png)

    ## saving PCA results to TSV...

    ## calculating euclidean distance between samples...

    ## generating sample distance heatmap...

    ## saving distance matrix to TSV...

![](deseq2pip_files/figure-html/unnamed-chunk-5-11.png)

    ## 
    ## ########################################################
    ## Running Differential Expression Analysis Pipeline
    ## ########################################################

    ## running differential expression analysis for the following comparisons:

    ## <<IgG_vs_IgM>>

    ##  - fitting DESeq2 model...

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 1303 DE genes found: 656 upregulated and 647 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-12.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-13.png)

    ## <<IgA_vs_IgM>>

    ##  - fitting DESeq2 model...

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 1723 DE genes found: 896 upregulated and 827 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-14.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-15.png)

    ## <<IgA_vs_IgG>>

    ##  - fitting DESeq2 model...

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 1279 DE genes found: 649 upregulated and 630 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-16.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-17.png)

    ## <<IgM_vs_IgG+IgA>>

    ##  - fitting DESeq2 model...

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 1049 DE genes found: 513 upregulated and 536 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-18.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-19.png)

    ## <<IgG_vs_IgM+IgA>>

    ##  - fitting DESeq2 model...

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 671 DE genes found: 339 upregulated and 332 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-20.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-21.png)

    ## <<IgA_vs_IgM+IgG>>

    ##  - fitting DESeq2 model...

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 1073 DE genes found: 528 upregulated and 545 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-22.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-5-23.png)

    ## 
    ## ########################################################
    ## Performing GSEA Analysis
    ## ########################################################

    ## <<IgG_vs_IgM>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-24.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-25.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-26.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-27.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-28.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-29.png)

    ## <<IgA_vs_IgM>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-30.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-31.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-32.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-33.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-34.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-35.png)

    ## <<IgA_vs_IgG>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-36.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-37.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-38.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-39.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-40.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-41.png)

    ## <<IgM_vs_IgG+IgA>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-42.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-43.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-44.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-45.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-46.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-47.png)

    ## <<IgG_vs_IgM+IgA>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-48.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-49.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-50.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-51.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-52.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-53.png)

    ## <<IgA_vs_IgM+IgG>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-54.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-55.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-56.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-57.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-58.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-5-59.png)

    ## exporting data for EnrichmentMap

    ## RNA-seq pipeline complete; returning DESeq2 object...

``` r
dds # return dds object
```

    ## class: DESeqDataSet 
    ## dim: 14178 11 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14178): 0610005C13Rik 0610009B22Rik ... Zzz3 a
    ## rowData names(1): gene
    ## colnames(11): BM_IgA_1 BM_IgA_2 ... BM_IgM_3 BM_IgM_4
    ## colData names(6): sample Group1 ... sizeFactor condition

## Specifying Custom Comparisons

By default, `comparisons = NULL` auto-generates all pairwise and
one-to-all comparisons from the factor levels of `var`. You can instead
pass a character vector of specific comparisons to run.

### Comparison string format

Comparisons follow the format `"numerator_vs_denominator"`. A positive
`log2FoldChange` means a gene is upregulated in the **numerator** group
relative to the **denominator** group.

``` r
# Run only selected pairwise comparisons
dds <- run_rna_pip(
    dds,
    org = "mouse",
    var = "Group2",
    design = "~ Group3 + Group2",
    comparisons = c("IgA_vs_IgM", "IgG_vs_IgM"),
    save_dir = save_dir)
```

### Combined-group comparisons

Use `+` to pool multiple groups on one side of the comparison:

``` r
# IgA and IgG pooled vs IgM
dds <- run_rna_pip(
    dds,
    org = "mouse",
    var = "Group2",
    design = "~ Group3 + Group2",
    comparisons = c("IgA+IgG_vs_IgM"),
    save_dir = save_dir)
```

### Mix pairwise and combined-group comparisons

``` r
dds <- run_rna_pip(
    dds,
    org = "mouse",
    var = "Group2",
    design = "~ Group3 + Group2",
    comparisons = c("IgA_vs_IgM", "IgG_vs_IgM", "IgA+IgG_vs_IgM"),
    save_dir = save_dir)
```

## Run deseq2pip on ATAC-seq data

Create a directory for results, currently set to tempdir()

``` r
save_dir <- tempdir()
dir.create(save_dir, recursive = TRUE)
```

    ## Warning in dir.create(save_dir, recursive = TRUE): '/tmp/RtmpuxUX6f' already
    ## exists

Load DESeq2 object processed by nf-core/atacseq

``` r
rdata <- system.file("data", "GSE224512.dds.RData", package = "deseq2pip")
annotatePeaks <- gzfile(system.file("data", "GSE224512.annotatePeaks.txt.gz", package = "deseq2pip"))
dds <- import_nfcore_atac(rdata = rdata, annotatePeaks = annotatePeaks)
```

Setup dds comparisons and design.

``` r
# subset groups for brevity
dds <- dds[, dds$Group1 %in% c("WT", 'BC', 'BCK')]

# set correct group factors
dds$Group1 <- factor(dds$Group1, c("WT", 'BC', 'BCK'))
dds$BC <- ifelse(dds$Group1 == "WT", "noBC", "BC")
dds$BC <- factor(dds$BC, c("noBC", "BC"))
dds$K <- ifelse(dds$Group1 == "BCK", "K", "noK")
dds$K <- factor(dds$K, c("noK", "K"))

# set design
design(dds) <- ~ BC + K
```

Run ATAC pipeline.

``` r
dds <- run_atac_pip(
    dds,
    org = "mouse",
    var = "BC",
    design = "~ BC + K",
    batch = NULL,
    comparisons = NULL,
    remove_xy = TRUE,
    remove_mt = TRUE,
    min_count = 10,
    order = "pxfc",
    TSS = TRUE,
    save_dir = save_dir)
```

    ## renv already initialized, restoring renv

    ## - The library is already synchronized with the lockfile.

    ## Running ATAC-seq pipeline with DESeq2pip v2.1.0

    ## 
    ## ########################################################
    ## Running Quality Control Pipeline
    ## ########################################################

    ## removing 3331 XY genes out of 74851 total genes...

    ## removing 2 MT genes out of 71520 total genes...

![](deseq2pip_files/figure-html/unnamed-chunk-12-1.png)

    ## after removing low expression genes, 68356 out of 71518 genes remain (95.58%)

![](deseq2pip_files/figure-html/unnamed-chunk-12-2.png)

    ## generating boxplots to check library size distribution...

![](deseq2pip_files/figure-html/unnamed-chunk-12-3.png)

    ## saving DESeq2 object & expressions...

    ## 
    ## ########################################################
    ## Running Distance Analysis Pipeline
    ## ########################################################

    ## performing variance stabilizing transformation (VST)...

    ## running PCA analysis...

    ## using ntop=500 top features by variance

    ##  - generating PCA plot for BC...

![](deseq2pip_files/figure-html/unnamed-chunk-12-4.png)

    ##  - generating PCA plot for Group1...

![](deseq2pip_files/figure-html/unnamed-chunk-12-5.png)

    ##  - generating PCA plot for Group2...

![](deseq2pip_files/figure-html/unnamed-chunk-12-6.png)

    ## saving PCA results to TSV...

    ## calculating euclidean distance between samples...

    ## generating sample distance heatmap...

    ## saving distance matrix to TSV...

![](deseq2pip_files/figure-html/unnamed-chunk-12-7.png)

    ## 
    ## ########################################################
    ## Running Differential Expression Analysis Pipeline
    ## ########################################################

    ## running differential expression analysis for the following comparisons:

    ## <<BC_vs_noBC>>

    ##  - fitting DESeq2 model...

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 299 DE genes found: 117 upregulated and 182 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-12-8.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-12-9.png)

    ##  - extracting annotations for all differentially expressed peaks...

    ##  - generating pie charts for 159 DE Peaks...

    ##  - generating pie charts for 62 Upregulated DE Peaks...

    ##  - generating pie charts for 97 Downregulated DE Peaks...

![](deseq2pip_files/figure-html/unnamed-chunk-12-10.png)![](deseq2pip_files/figure-html/unnamed-chunk-12-11.png)![](deseq2pip_files/figure-html/unnamed-chunk-12-12.png)

    ## 
    ## ########################################################
    ## Running Differential Expression Analysis Pipeline
    ## ########################################################

    ## running differential expression analysis for the following comparisons:

    ## <<BC_vs_noBC>>

    ##  - fitting DESeq2 model...

    ##  - extracting differential expression results...

    ##  - shrinking lfc by ashr...

    ##  - 62 DE genes found: 28 upregulated and 34 downregulated

    ##  - generating MA plot...

![](deseq2pip_files/figure-html/unnamed-chunk-12-13.png)

    ##  - generating volcano plot...

![](deseq2pip_files/figure-html/unnamed-chunk-12-14.png)

    ## 
    ## ########################################################
    ## Performing GSEA Analysis
    ## ########################################################

    ## <<BC_vs_noBC>>

    ##  - running GSEA for BIOCARTA...

    ##  - saving BIOCARTA results to RDS and TSV...

    ##  - running GSEA for GOBP...

    ##  - saving GOBP results to RDS and TSV...

    ##  - running GSEA for HALLMARK...

    ##  - saving HALLMARK results to RDS and TSV...

    ##  - running GSEA for KEGG...

    ##  - saving KEGG results to RDS and TSV...

    ##  - running GSEA for REACTOME...

    ##  - saving REACTOME results to RDS and TSV...

    ##  - running GSEA for TFT...

    ##  - saving TFT results to RDS and TSV...

    ##  - generating top 10 BIOCARTA barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-12-15.png)

    ##  - generating top 10 GOBP barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-12-16.png)

    ##  - generating top 10 HALLMARK barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-12-17.png)

    ##  - generating top 10 KEGG barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-12-18.png)

    ##  - generating top 10 REACTOME barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-12-19.png)

    ##  - generating top 10 TFT barplots...

![](deseq2pip_files/figure-html/unnamed-chunk-12-20.png)

    ## exporting data for EnrichmentMap

    ## ATAC-seq pipeline complete; returning DESeq2 object...

``` r
dds # return dds object
```

    ## class: DESeqDataSet 
    ## dim: 68356 13 
    ## metadata(1): version
    ## assays(2): counts vst
    ## rownames(68356): Interval_1 Interval_2 ... Interval_74799
    ##   Interval_74821
    ## rowData names(20): peaks Chr ... gene TSS
    ## colnames(13): WT_REP4 BCK_REP5 ... BC_REP2 BC_REP1
    ## colData names(6): sample Group1 ... BC K
