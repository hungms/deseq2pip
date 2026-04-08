# Read Differential Expression Results

This function reads differential expression results from multiple
comparison directories and optionally merges them into a single data
frame.

This function reads differential expression results from multiple
comparison directories and optionally merges them into a single data
frame.

## Usage

``` r
read_diffexp(group_dir = getwd(), merge = T)

read_diffexp(group_dir = getwd(), merge = T)
```

## Arguments

- group_dir:

  Parent directory containing comparison subdirectories. Default is
  current working directory

- merge:

  Logical. If TRUE, returns a single merged data frame. If FALSE,
  returns a list of data frames

## Value

Either a merged data frame or a list of data frames containing
differential expression results

Either a merged data frame or a list of data frames containing
differential expression results

## Examples

``` r
if (FALSE) { # \dontrun{
# Read and merge all results
merged_results <- read_diffexp(group_dir = "results", merge = TRUE)

# Read results as a list
results_list <- read_diffexp(group_dir = "results", merge = FALSE)
} # }
if (FALSE) { # \dontrun{
# Read and merge all results
merged_results <- read_diffexp(group_dir = "results", merge = TRUE)

# Read results as a list
results_list <- read_diffexp(group_dir = "results", merge = FALSE)
} # }
```
