# Read GSEA Results from RDS Files

This function reads GSEA results from RDS files for multiple comparisons
and collections.

## Usage

``` r
read_gsea_rds(
  group_dir = getwd(),
  collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT")
)
```

## Arguments

- group_dir:

  Parent directory containing comparison subdirectories. Default is
  current working directory

- collection:

  Vector of gene set collections to load. Default includes HALLMARK,
  GOBP, KEGG, REACTOME, BIOCARTA, and TFT

## Value

A list of GSEA objects for each comparison

## Examples

``` r
if (FALSE) { # \dontrun{
# Read all default collections
gsea_results <- read_gsea_rds("results")

# Read specific collections
gsea_results <- read_gsea_rds("results", collection = c("HALLMARK", "KEGG"))
} # }
```
