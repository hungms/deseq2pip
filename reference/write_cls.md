# Write class file pipeline for EnrichmentMap

This function writes a class file for EnrichmentMap visualization.

## Usage

``` r
write_cls(dds, var, cls_name = "dds_class.cls", save_dir = getwd())
```

## Arguments

- dds:

  DESeq2 object

- var:

  Group by column name

- cls_name:

  Name of the class file. Default is "dds_class.cls"

- save_dir:

  Directory where the formatted files will be saved. Default is current
  working directory
