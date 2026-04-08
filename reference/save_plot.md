# Save ggplot Object as PDF

This function saves a ggplot object as a PDF file with specified
dimensions. It tries to use Cairo PDF device if available, or falls back
to standard PDF if not.

## Usage

``` r
save_plot(input, plot_name, save_dir, w, h)
```

## Arguments

- input:

  ggplot object to be saved

- plot_name:

  Name of the output PDF file

- save_dir:

  Directory where the PDF file will be saved

- w:

  Width of the PDF in inches

- h:

  Height of the PDF in inches

## Value

The file path (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
# Save plot with default dimensions
save_plot(my_plot, "plot.pdf", "figures", w = 8, h = 6)

# Save plot with custom dimensions
save_plot(my_plot, "plot.pdf", "figures", w = 10, h = 8)
} # }
```
