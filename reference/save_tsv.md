# Save Data Frame as TSV File

This function saves a data frame as a tab-separated values (TSV) file.

## Usage

``` r
save_tsv(input, tsv_name, save_dir, row.names = F)
```

## Arguments

- input:

  Data frame to be saved

- tsv_name:

  Name of the output TSV file

- save_dir:

  Directory where the TSV file will be saved

- row.names:

  Logical. If TRUE, row names will be included in the output. Default is
  FALSE

## Value

None. Creates a TSV file in the specified directory

## Examples

``` r
if (FALSE) { # \dontrun{
# Save data frame without row names
save_tsv(my_data, "output.tsv", "results")

# Save data frame with row names
save_tsv(my_data, "output.tsv", "results", row.names = TRUE)
} # }
```
