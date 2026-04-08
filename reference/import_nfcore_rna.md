# Import nfcore/rnaseq DESeq2 Object

This function imports a DESeq2 object from nfcore RNA-seq pipeline
output and adds gene names.

## Usage

``` r
import_nfcore_rna(rdata, tx2gene)
```

## Arguments

- rdata:

  Path to the RData file containing the DESeq2 object

- tx2gene:

  Path to the tx2gene mapping file

## Value

A DESeq2 object with gene names added to rowData

## Examples

``` r
if (FALSE) { # \dontrun{
# Import DESeq2 object and add gene names
dds <- import_nfcore_rna("path/to/dds.RData", "path/to/tx2gene.tsv")
} # }
```
