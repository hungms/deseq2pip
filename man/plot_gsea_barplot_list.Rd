% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea_functions.R
\name{plot_gsea_barplot_list}
\alias{plot_gsea_barplot_list}
\title{Generate GSEA Barplots for All Comparisons}
\usage{
plot_gsea_barplot_list(gsea.df, ...)
}
\arguments{
\item{gsea.df}{GSEA result data frame from run_gsea()}

\item{...}{Additional arguments passed to plot_gsea_barplot()}
}
\value{
A list of GSEA barplots for each comparison and collection
}
\description{
This function creates barplots to visualize GSEA results for all comparisons
and gene set collections.
}
\examples{
\dontrun{
# Generate GSEA barplots for all results
plots <- plot_gsea_barplot_list(gsea_results)

# Generate with custom parameters
plots <- plot_gsea_barplot_list(gsea_results, n = 15, signif = TRUE)
}
}
