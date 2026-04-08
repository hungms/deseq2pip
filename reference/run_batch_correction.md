# Run batch correction by Limma & Combat

This function performs batch correction on the variance-stabilized
transformed data using the Limma package. It can help identify potential
batch effects or outliers.

## Usage

``` r
run_batch_correction(
  dds,
  batch,
  method = c("limma", "combat"),
  save_data = TRUE,
  save_dir = getwd()
)
```

## Arguments

- dds:

  A DESeq2 object containing the normalized gene expression data to be
  batch corrected

- batch:

  Column name in colData(vsd) to use for batch correction.

- method:

  Method to use for batch correction. Options are "limma" and "combat".

- save_data:

  Logical. If TRUE, saves the batch-corrected data to a TSV file.
  Default is TRUE

- save_dir:

  Directory to save the results. Default is the current working
  directory

## Value

A DESeq2 object containing the batch-corrected data
