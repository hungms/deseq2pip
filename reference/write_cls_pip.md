# Write class file pipeline for EnrichmentMap

This function writes a class file for EnrichmentMap visualization.

## Usage

``` r
write_cls_pip(dds, var, group_dir = group_dir, save_dir = getwd())
```

## Arguments

- dds:

  DESeq2 object

- var:

  Group by column name

- group_dir:

  Parent directory containing comparison subdirectories. Default is
  current working directory

- save_dir:

  Directory where the formatted files will be saved. Default is current
  working directory
