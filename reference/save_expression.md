# Save Expression Data from DESeq2 Object

This function saves expression data from a DESeq2 object in various
formats, including the DESeq2 object itself, raw counts, normalized
expression values, and class labels for GSEA.

## Usage

``` r
save_expression(dds, var, save_dir = getwd())
```

## Arguments

- dds:

  DESeq2 object containing the expression data

- save_dir:

  Directory to save files. Default is the current working directory

- save_dir_name:

  Name of the subdirectory to save files in. Default is "qc_results"

## Value

Nothing, but saves the following files: - DESeq2 object as RDS - Raw
counts as TSV - Variance-stabilized transformed data as TSV - CLS file
with class labels

## Examples

``` r
if (FALSE) { # \dontrun{
# Save expression data with default settings
save_expression(dds)

# Save expression data with custom grouping and directory name
save_expression(dds, var = "Treatment", save_dir_name = "custom_results")
} # }
```
