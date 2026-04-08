# Plot ATAC-seq Peak Annotation

This function plots the annotation distribution of ATAC-seq peaks.

## Usage

``` r
plot_peak_annot_pip(dds, res, save_dir = getwd(), ...)
```

## Arguments

- dds:

  DESeq2 object containing the ATAC-seq peaks data

- res:

  Differential expression results data frame from run_diffexp()

- save_dir:

  Directory to save the plot. Default is current working directory

- ...:

  Additional arguments (ignored, for compatibility)

## Value

A list of ggplot objects showing the pie charts
