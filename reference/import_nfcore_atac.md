# Import nfcore ATAC-seq DESeq2 Object

This function imports a DESeq2 object from nfcore ATAC-seq pipeline
output and adds gene names.

## Usage

``` r
import_nfcore_atac(rdata, annotatePeaks, dist.to.TSS = 2000)
```

## Arguments

- rdata:

  Path to the RData file containing the DESeq2 object

- dist.to.TSS:

  Distance to TSS to consider as TSS. Default is 2000

- tx2gene:

  Path to the tx2gene mapping file

## Value

A DESeq2 object with gene names added to rowData

## Examples

``` r
if (FALSE) { # \dontrun{
# Import DESeq2 object and add gene names
dds <- import_nfcore_atac("path/to/dds.RData", "path/to/tx2gene.tsv")
} # }
```
