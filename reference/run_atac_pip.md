# Run Complete ATAC-seq Pipeline

This function runs a complete ATAC-seq analysis pipeline including
quality control, differential expression analysis, and gene set
enrichment analysis for all possible comparisons in the experiment. It
generates various plots and saves results in organized directories.

## Usage

``` r
run_atac_pip(
  dds,
  org,
  var,
  design,
  batch = NULL,
  comparisons = NULL,
  pals = NULL,
  remove_xy = FALSE,
  remove_mt = FALSE,
  min_count = 10,
  order = "pxfc",
  TSS = TRUE,
  save_dir = getwd()
)
```

## Arguments

- dds:

  DESeq2 object

- org:

  The organism to use, either "human" or "mouse".

- var:

  Column name in colData(dds) to use for grouping.

- design:

  Design formula for the DESeq2 object.

- comparisons:

  Vector of comparisons to run. Default is NULL

- pals:

  Vector of colors to use for groups. If NULL, uses default ggplot2
  colors. Default is NULL

- remove_xy:

  Logical. If TRUE, removes genes on X and Y chromosomes. Default is
  FALSE

- remove_mt:

  Logical. If TRUE, removes mitochondrial genes. Default is FALSE

- min_count:

  Minimum count threshold for filtering lowly expressed genes. Default
  is 10

- order:

  Column name to use for ranking genes. Default is "pxfc"

- TSS:

  Logical. If TRUE, repeats the analysis for TSS peaks. Default is TRUE

- save_dir:

  Directory where all output files will be saved. Default is current
  working directory

## Value

The processed DESeq2 object

## Examples

``` r
if (FALSE) { # \dontrun{
# Run complete pipeline with default settings
dds <- run_atac_pip(dds)

# Run pipeline with custom settings
dds <- run_atac_pip(dds, 
                    org = "mouse", 
                    remove_xy = TRUE,
                    var = "Treatment",
                    save_dir = "output")

} # }
```
