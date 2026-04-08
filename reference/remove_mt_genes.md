# Remove Mitochondrial Genes from DESeq2 Object

This function removes genes encoded in the mitochondrial genome from a
DESeq2 object. It uses the Ensembl database to identify mitochondrial
genes.

## Usage

``` r
remove_mt_genes(dds, org, ...)
```

## Arguments

- dds:

  A DESeq2 object containing the gene expression data

- org:

  The organism to use, either "human" or "mouse".

## Value

A filtered DESeq2 object with mitochondrial genes removed

## Examples

``` r
if (FALSE) { # \dontrun{
# For human data
dds_filtered <- remove_mt_genes(dds, org)

# For mouse data
dds_filtered <- remove_mt_genes(dds, org = "mouse")
} # }
```
