# Remove Genes with Low Expression

This function filters out genes with raw count values below a specified
minimum count threshold in at least a minimum number of replicates per
condition.

## Usage

``` r
remove_low_expression(
  dds,
  var,
  min_count = 10,
  save_plot = TRUE,
  save_dir = getwd()
)
```

## Arguments

- dds:

  A DESeq2 object containing the gene expression data

- var:

  Column name in colData(dds) to use for defining conditions. Default is
  "Group1"

- min_count:

  Minimum count threshold for filtering. Default is 10

- save_plot:

  Logical. If TRUE, saves the expression distribution plot. Default is
  TRUE

- save_dir:

  Directory to save the plot. Default is the current working directory

## Value

A filtered DESeq2 object with low-expression genes removed

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove genes with counts < 10
dds_filtered <- remove_low_expression(dds, min_count = 10)

# Remove genes with counts < 5 and save plot
dds_filtered <- remove_low_expression(dds, min_count = 5, save_plot = TRUE)
} # }
```
