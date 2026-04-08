# Log output

This function logs the script to the specified directory.

## Usage

``` r
log_output(call, save_dir, expr = NULL)
```

## Arguments

- save_dir:

  Directory to save the script. Default is the current working directory

- expr:

  Expression to evaluate and log. If NULL, only logs the function call
  and arguments

## Value

Nothing, but logs the script to the specified directory

## Examples

``` r
if (FALSE) { # \dontrun{
# Log script
log_script(save_dir = "results")
} # }
```
