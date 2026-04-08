# Format Enrichment Map Data

This function processes and formats GSEA results for visualization in
EnrichmentMap.

## Usage

``` r
enrichmentmap_pip(
  dds,
  org,
  var,
  collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME"),
  group_dir,
  save_dir = group_dir,
  save_dir_name = "enrichmentmap"
)
```

## Arguments

- dds:

  DESeq2 object

- org:

  The organism to use, either "human" or "mouse".

- collection:

  Vector of gene set collections to process. Default includes HALLMARK,
  GOBP, KEGG, and REACTOME

- group_dir:

  Parent directory containing comparison subdirectories. Default is
  current working directory

- save_dir:

  Directory where the formatted files will be saved. Default is current
  working directory

- save_dir_name:

  Name of the subdirectory to save files in. Default is "enrichmentmap"

## Value

None. Creates formatted files for EnrichmentMap visualization

## Examples

``` r
if (FALSE) { # \dontrun{
# Format all default collections
enrichmentmap_pip(dds, var = "group", group_dir = "results")

# Format specific collections
enrichmentmap_pip(dds, var = "group", group_dir = "results", collection = c("HALLMARK", "KEGG"))
} # }
```
