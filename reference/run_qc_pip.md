# Run Quality Control Pipeline

This function runs the quality control portion of the RNA-seq analysis
pipeline, including filtering genes, generating QC plots, and saving
expression data.

## Usage

``` r
run_qc_pip(
  dds,
  org,
  var,
  remove_xy = TRUE,
  remove_mt = TRUE,
  min_count = 10,
  save_dir = getwd(),
  save_dir_name = "qc_results"
)
```

## Arguments

- dds:

  DESeq2 object containing the expression data

- org:

  The organism to use, either "human" or "mouse".

- var:

  Column name in colData(dds) to use for grouping.

- remove_xy:

  Logical. If TRUE, removes genes on X and Y chromosomes. Default is
  TRUE

- remove_mt:

  Logical. If TRUE, removes mitochondrial genes. Default is TRUE

- min_count:

  Minimum count threshold for filtering lowly expressed genes. Default
  is 10

- save_dir:

  Directory where all output files will be saved. Default is current
  working directory

- save_dir_name:

  Name of the subdirectory to save files in. Default is "qc_results"

## Value

The processed DESeq2 object

## Examples

``` r
if (FALSE) { # \dontrun{
# Run QC pipeline with default settings
dds <- run_qc_pip(dds)

# Run QC pipeline with custom settings
dds <- run_qc_pip(dds, org = "mouse", remove_xy = TRUE, save_dir_name = "custom_results")
} # }
```
