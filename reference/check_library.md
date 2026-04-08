# Check Library Size Distribution

This function generates a boxplot showing the distribution of expression
values across samples. It helps identify potential outliers or quality
issues in the sequencing libraries.

## Usage

``` r
check_library(dds, save_plot = TRUE, save_dir = getwd())
```

## Arguments

- dds:

  A DESeq2 object containing the gene expression data

- save_plot:

  Logical. If TRUE, saves the library size plot to PDF. Default is TRUE

- save_dir:

  Directory to save the plot. Default is the current working directory

## Value

A ggplot object showing the library size distribution

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate and display plot
p <- check_library(dds)

# Generate and save plot
p <- check_library(dds, save_plot = TRUE, save_dir_name = "custom_results")
} # }
```
