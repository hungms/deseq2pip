#' Package Load Function
#'
#' This function is called when the package is loaded. It performs the following tasks:
#' 1. Loads required packages for the DESeq2 pipeline
#' 2. Sets up default options for dplyr and random number generation
#' 3. Initializes error message handling
#'
#' @param libname The library directory where the package is installed
#' @param pkgname The name of the package
#' @return None. This function is called for its side effects
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
    "strpip", "plotr", "DESeq2", "biomaRt", "ggalt", "pheatmap", "ggplotify", "clusterProfiler", "msigdbr", "fgsea", "ggpubr", "grid", "scales", "forcats")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
  errorMessage <- NULL

  #for(o in c("human", "mouse")){
  #  assign(paste0(o, "_msigdb"), import_msigdb(org = o), envir = parent.env(environment()))}

  ### load essential scripts
  #dir <- system.file("essentials", package = pkgname)
  #scripts <- list.files(dir, full.names = TRUE)
  #for(script in scripts){
  #  source(script)}

}
