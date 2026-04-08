# Get TSS peaks from DESeq2 object

This function extracts TSS peaks from a DESeq2 object and removes
duplicate genes.

## Usage

``` r
getTSS(dds)
```

## Arguments

- dds:

  A DESeq2 object

## Value

A DESeq2 object containing only TSS peaks

## Examples

``` r
if (FALSE) { # \dontrun{
# Get TSS peaks
dds.tss <- getTSS(dds)
} # }
```
