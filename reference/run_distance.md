# Calculate and Plot Sample Distances

This function calculates the Euclidean distance between samples based on
their expression profiles and generates a heatmap visualization. It
helps identify sample relationships and potential outliers.

## Usage

``` r
run_distance(
  vsd,
  save_prefix = NULL,
  save_data = TRUE,
  save_plot = TRUE,
  save_dir = getwd(),
  ...
)
```

## Arguments

- vsd:

  A DESeq2 object containing the normalized gene expression data

- save_prefix:

  Prefix for the save file. Default is NULL

- save_data:

  Logical. If TRUE, saves distance matrix to TSV. Default is TRUE

- save_plot:

  Logical. If TRUE, saves the distance heatmap to PDF. Default is TRUE

- save_dir:

  Directory to save the results. Default is the current working
  directory

- ...:

  Additional arguments passed to pheatmap()

## Value

A pheatmap object showing the sample distances

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate and display heatmap
p <- run_distance(vsd)

# Generate and save heatmap with custom parameters
p <- run_distance(vsd, save_plot = TRUE, save_dir_name = "custom_results")
} # }
```
