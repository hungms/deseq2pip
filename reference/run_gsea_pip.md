# Run GSEA pipeline

This function runs the GSEA pipeline for all comparisons and
collections.

## Usage

``` r
run_gsea_pip(
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

- msigdbr:

  Dataframe or path to a custom gene set database in GMT format. Default
  is NULL

- save_data:

  Logical. If TRUE, saves results to RDS and TSV files. Default is TRUE

- group_save_dir:

  Directory to save the results.

## Value

A list of GSEA results for each gene set collection
