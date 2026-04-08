# Plot GSEA Results

This function creates a barplot to visualize GSEA results, showing the
most significant gene sets in each direction (up and down-regulated).

## Usage

``` r
plot_gsea_barplot(gsea, n = 10, signif = F, save_plot = T, save_dir = getwd())
```

## Arguments

- gsea:

  GSEA result data frame from run_gsea()

- n:

  Number of top gene sets to show in each direction. Default is 10

- signif:

  Logical. If TRUE, only shows gene sets with q-value \< 0.05. Default
  is TRUE

- save_plot:

  Logical. If TRUE, saves the plot to PDF. Default is TRUE

- save_dir:

  Directory to save the plot. Default is current working directory

## Value

A ggplot object showing the GSEA barplot

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic GSEA plot
p <- plot_gsea_barplot(gsea_results)

# GSEA plot with more gene sets and only significant results
p <- plot_gsea_barplot(gsea_results, n = 15, signif = TRUE)
} # }
```
