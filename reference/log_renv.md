# Log renv Snapshot

This function logs the renv snapshot to the specified directory.

## Usage

``` r
log_renv(save_dir)
```

## Arguments

- save_dir:

  Directory to save the renv snapshot. Default is the current working
  directory

## Value

Nothing, but logs the renv snapshot to the specified directory

## Examples

``` r
if (FALSE) { # \dontrun{
# Log renv snapshot
log_renv(save_dir = "results")
} # }
```
