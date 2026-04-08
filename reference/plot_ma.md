# Generate MA Plot

This function creates a MA plot to visualize differential expression
results. It shows the relationship between mean expression and log2 fold
change for each gene.

## Usage

``` r
plot_ma(res, fc.thresh = 0.5, save_plot = TRUE, save_dir = getwd())
```

## Arguments

- res:

  Differential expression result data frame from run_diffexp()

- save_plot:

  Logical. If TRUE, saves the plot to PDF. Default is TRUE

- save_dir:

  Directory to save the plot. Default is current working directory

## Value

A ggplot object showing the MA plot

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic MA plot
p <- plot_ma(res)
} # }
```
