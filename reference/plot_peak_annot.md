# Plot ATAC-seq Peak Annotations

This function analyzes differentially expressed ATAC-seq peaks and
creates pie charts of their genomic annotations. It filters peaks by
adjusted p-value and log2 fold change, cleans up annotation labels, and
generates three pie charts: one for all DE peaks, one for upregulated
peaks, and one for downregulated peaks.

## Usage

``` r
plot_peak_annot(
  dds,
  res,
  p.thresh = 0.05,
  fc.thresh = 0.5,
  save_plot = TRUE,
  save_dir = getwd()
)
```

## Arguments

- dds:

  DESeq2 object containing the ATAC-seq peaks data

- res:

  Differential expression results data frame from run_diffexp()

- p.thresh:

  Adjusted p-value threshold for significance. Default is 0.05

- fc.thresh:

  Log2 fold change threshold for significance (absolute value). Default
  is 0.5

- save_plot:

  Logical. If TRUE, saves the plots to PDF. Default is TRUE

- save_dir:

  Directory to save the plot. Default is current working directory

## Value

A list of ggplot objects showing the pie charts

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate annotation pie charts for ATAC-seq peaks
plots <- plot_peak_annot(dds, res)

# Use custom thresholds
plots <- plot_peak_annot(dds, res, p.thresh = 0.01, fc.thresh = 1)
} # }
```
