# Generate Volcano Plot

This function creates a volcano plot to visualize differential
expression results. It shows the relationship between statistical
significance (-log10 adjusted p-value) and biological significance (log2
fold change) for each gene.

## Usage

``` r
plot_volcano(
  res,
  n = 25,
  fc.thresh = 0.5,
  p.thresh = 0.05,
  crop = T,
  highlight.genes = NULL,
  save_plot = T,
  save_dir = getwd()
)
```

## Arguments

- res:

  Differential expression result data frame from run_diffexp()

- n:

  Number of top genes to label in each direction. Default is 25

- fc.thresh:

  Log2 fold change threshold for significance. Default is 1

- p.thresh:

  Adjusted p-value threshold for significance. Default is 0.05

- crop:

  Logical. If TRUE, limits the x-axis range to genes with significant
  changes. Default is TRUE

- highlight.genes:

  Vector of gene names to highlight. Default is NULL

- save_plot:

  Logical. If TRUE, saves the plot to PDF. Default is TRUE

- save_dir:

  Directory to save the plot. Default is current working directory

## Value

A ggplot object showing the volcano plot

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic volcano plot
p <- plot_volcano(res)

# Volcano plot with custom thresholds and highlighted genes
p <- plot_volcano(res, n = 30, fc.thresh = 2, p.thresh = 0.01,
                  highlight.genes = c("GENE1", "GENE2"))
} # }
```
