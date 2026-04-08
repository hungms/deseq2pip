# Run Differential Expression Wrapper

This function runs the differential expression analysis and generates MA
and volcano plots for a single comparison

## Usage

``` r
run_diffexp_pip(
  dds,
  design,
  org,
  var,
  comparisons = NULL,
  order = "pxfc",
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

- comparisons:

  Comparison to run.

- order:

  Column name to use for ranking genes. Default is "pxfc"

- one_to_all:

  Logical. If TRUE, runs one-to-all comparisons. Default is FALSE

- group_save_dir:

  Directory where all output files will be saved. Default is current
  working directory

## Value

A data frame containing differential expression results for a comparison

## Examples

``` r
if (FALSE) { # \dontrun{
# Run differential expression pipeline
res <- run_diffexp_wrapper(dds)

# Run with custom grouping
res <- run_diffexp_wrapper(dds, var = "Treatment")
} # }
```
