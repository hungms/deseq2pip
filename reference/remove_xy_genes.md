# Remove XY Chromosome Genes from DESeq2 Object

This function removes genes located on the X and Y chromosomes from a
DESeq2 object. It uses the Ensembl database to identify genes on these
chromosomes.

## Usage

``` r
remove_xy_genes(dds, org, ...)
```

## Arguments

- dds:

  A DESeq2 object containing the gene expression data

- org:

  The organism to use, either "human" or "mouse".

## Value

A filtered DESeq2 object with XY chromosome genes removed

## Examples

``` r
if (FALSE) { # \dontrun{
# For human data
dds_filtered <- remove_xy_genes(dds, org)

# For mouse data
dds_filtered <- remove_xy_genes(dds, org = "mouse")
} # }
```
