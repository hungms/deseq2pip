# Import MSigDB Gene Sets

This function imports pre-defined MSigDB gene sets for human or mouse
organisms.

## Usage

``` r
import_msigdbr(org)
```

## Arguments

- org:

  The organism to use, either "human" or "mouse".

## Value

A data frame containing MSigDB gene sets

## Examples

``` r
if (FALSE) { # \dontrun{
# Import human MSigDB gene sets
msigdbr_human <- import_msigdbr("human")

# Import mouse MSigDB gene sets
msigdbr_mouse <- import_msigdbr("mouse")
} # }
```
