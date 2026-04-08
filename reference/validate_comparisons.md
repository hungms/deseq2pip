# Validate comparisons

This function validates that all comparisons listed are found in
dds\[[var](https://rdrr.io/r/stats/cor.html)\]. If comparisons is NULL,
it will automatically generate all possible comparisons.

## Usage

``` r
validate_comparisons(dds, var, comparisons)
```

## Arguments

- dds:

  A DESeq2 object.

- var:

  Column name in colData(dds) to use for grouping.

- comparisons:

  Vector of comparisons to validate. If NULL, generates all possible
  comparisons.

## Value

A validated vector of comparisons.
