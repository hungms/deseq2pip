# Validate design formula

This function validates a design formula and ensures all colData columns
used in the design are factors with levels set.

## Usage

``` r
validate_design(dds, design)
```

## Arguments

- dds:

  A DESeq2 object.

- design:

  A design formula (character string or formula object).

## Value

A DESeq2 object with validated and converted factors.
