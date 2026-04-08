# Run Principal Component Analysis

This function performs PCA on the variance-stabilized transformed data
and generates a PCA plot. It can visualize sample relationships and
identify potential batch effects or outliers.

## Usage

``` r
run_pca(
  vsd,
  var,
  pals = NULL,
  size = 4,
  save_prefix = NULL,
  save_data = TRUE,
  save_plot = TRUE,
  save_dir = getwd()
)
```

## Arguments

- vsd:

  A DESeq2 object containing the normalized gene expression data

- var:

  Column names in colData(vsd) to group by.

- pals:

  Vector of colors to use for groups. If NULL, uses default ggplot2
  colors. Default is NULL

- size:

  Size of points in the PCA plot. Default is 4

- save_prefix:

  Prefix for the save file. Default is NULL

- save_data:

  Logical. If TRUE, saves PCA results to TSV. Default is TRUE

- save_plot:

  Logical. If TRUE, saves the PCA plot to PDF. Default is TRUE

- save_dir:

  Directory to save the results. Default is the current working
  directory

- shape:

  Column name in colData(vsd) to use for shape in the PCA plot. Default
  is NULL

## Value

A ggplot object showing the PCA plot

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic PCA plot
p <- run_pca(vsd)

# PCA plot with shape and save results
p <- run_pca(vsd, shape = "Batch", save_data = TRUE)

# PCA plot with custom colors
p <- run_pca(vsd, pals = c("red", "blue", "green"))
} # }
```
