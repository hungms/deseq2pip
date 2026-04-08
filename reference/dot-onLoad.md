# Package Load Function

This function is called when the package is loaded. It performs the
following tasks:

1.  Loads required packages for the DESeq2 pipeline

2.  Sets up default options for dplyr and random number generation

3.  Initializes error message handling

## Usage

``` r
.onLoad(libname, pkgname)
```

## Arguments

- libname:

  The library directory where the package is installed

- pkgname:

  The name of the package

## Value

None. This function is called for its side effects
