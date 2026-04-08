# Read GSEA Results from TSV Files

This function reads GSEA results from TSV files for multiple comparisons
and collections.

## Usage

``` r
read_gsea_tsv(
  group_dir = getwd(),
  collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT"),
  merge = T
)
```

## Arguments

- group_dir:

  Parent directory containing comparison subdirectories. Default is
  current working directory

- collection:

  Vector of gene set collections to load. Default includes HALLMARK,
  GOBP, KEGG, REACTOME, BIOCARTA, and TFT

- merge:

  Logical. If TRUE, returns a single merged data frame. If FALSE,
  returns a list of data frames

## Value

Either a merged data frame or a list of data frames containing GSEA
results

## Examples

``` r
if (FALSE) { # \dontrun{
# Read and merge all results
merged_gsea <- read_gsea_tsv("results", merge = TRUE)

# Read specific collections as a list
gsea_list <- read_gsea_tsv("results", collection = c("HALLMARK", "KEGG"), merge = FALSE)
} # }
```
