# Run Differential Expression Analysis

This function runs the differential expression analysis for a single
comparison.

## Usage

``` r
run_diffexp(
  dds,
  org,
  var,
  design,
  comparison,
  order,
  save_data = TRUE,
  save_dir = getwd()
)
```

## Arguments

- dds:

  DESeq2 object containing the expression data

- org:

  The organism to use, either "human" or "mouse".

- var:

  Column name in colData(dds) to use for grouping.

- design:

  Design formula for the DESeq2 object.

- comparison:

  Comparison to run.

- order:

  Column name to use for ranking genes. Default is "pxfc"

- save_data:

  Logical. If TRUE, saves the results to a TSV file. Default is TRUE

- save_dir:

  Directory to save the results. Default is current working directory

## Value

A data frame containing differential expression results for a comparison
