# Run Sample Distance Pipeline This function runs the distance analysis portion of the RNA-seq analysis pipeline, including generating PCA and distance plots, and saving distance data.

Run Sample Distance Pipeline This function runs the distance analysis
portion of the RNA-seq analysis pipeline, including generating PCA and
distance plots, and saving distance data.

## Usage

``` r
run_dist_pip(
  dds,
  var,
  batch = NULL,
  pals = NULL,
  save_dir = getwd(),
  save_dir_name = "qc_results"
)
```

## Arguments

- dds:

  DESeq2 object containing the expression data

- var:

  Column name in colData(dds) to use for grouping. Default is "Group1"

- pals:

  Vector of colors to use for groups. If NULL, uses default ggplot2
  colors. Default is NULL

- save_dir:

  Directory where all output files will be saved. Default is current
  working directory

- save_dir_name:

  Name of the subdirectory to save files in. Default is "qc_results"

- design:

  Design formula for the DESeq2 object.

## Value

The processed DESeq2 object

## Examples

``` r
if (FALSE) { # \dontrun{
# Run distance pipeline with default settings
dds <- run_dist_pip(dds)

# Run distance pipeline with custom settings and batch correction
dds <- run_dist_pip(dds, batch = "Batch", save_dir_name = "qc_results")
} # }
```
