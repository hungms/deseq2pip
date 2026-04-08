# Run Gene Set Enrichment Analysis

This function performs Gene Set Enrichment Analysis (GSEA) using
clusterProfiler to identify enriched gene sets in the differential
expression results. It can use either pre-defined MSigDB gene sets or
custom gene sets.

## Usage

``` r
run_gsea(
  res,
  org,
  msigdbr = import_msigdbr(org),
  save_data = T,
  save_dir = getwd()
)
```

## Arguments

- res:

  Differential expression result data frame from run_diffexp()

- org:

  The organism to use, either "human" or "mouse".

- save_data:

  Logical. If TRUE, saves results to TSV and RDS files. Default is TRUE

- save_dir:

  Directory to save the results. Default is current working directory

- custom_msigdb:

  Dataframe or path to a custom gene set database in GMT format. Default
  is NULL

## Value

A list of GSEA results for each gene set collection, with each element
containing: - gene set name - normalized enrichment score (NES) -
p-value and adjusted p-value - leading edge genes - collection name -
comparison name

## Examples

``` r
if (FALSE) { # \dontrun{
# Run GSEA with default MSigDB gene sets
gsea_results <- run_gsea(res)

# Run GSEA with custom gene sets
gsea_results <- run_gsea(res, custom_msigdb = "path/to/custom_sets.tsv")

# Run GSEA for mouse data
gsea_results <- run_gsea(res, org = "mouse")
} # }
```
