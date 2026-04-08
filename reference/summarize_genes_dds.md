# Summarize Gene Expression in DESeq2 Object

This function aggregates expression values from multiple isoforms/peaks
of the same gene in a DESeq2 object into a single gene-level expression
value. It creates a new DESeq2 object with gene-level expression data.

## Usage

``` r
summarize_genes_dds(dds, gene_sym_col = "gene", ...)
```

## Arguments

- dds:

  DESeq2 object containing isoform/peak-level expression data

- gene_sym_col:

  Column name in rowData(dds) containing gene symbols. Default is
  "Gene.Name"

- ...:

  Additional arguments passed to summarize_genes()

## Value

A new DESeq2 object containing gene-level expression data

## Examples

``` r
if (FALSE) { # \dontrun{
# Sum gene expression
gene_dds <- summarize_genes_dds(dds)

# Average gene expression
gene_dds <- summarize_genes_dds(dds, normalized = TRUE)
} # }
```
