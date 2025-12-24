#' Save Data Frame as TSV File
#'
#' This function saves a data frame as a tab-separated values (TSV) file.
#'
#' @param input Data frame to be saved
#' @param tsv_name Name of the output TSV file
#' @param save_dir Directory where the TSV file will be saved
#' @param row.names Logical. If TRUE, row names will be included in the output. Default is FALSE
#' @return None. Creates a TSV file in the specified directory
#' @examples
#' \dontrun{
#' # Save data frame without row names
#' save_tsv(my_data, "output.tsv", "results")
#' 
#' # Save data frame with row names
#' save_tsv(my_data, "output.tsv", "results", row.names = TRUE)
#' }
#' @export
save_tsv <- function(input, tsv_name, save_dir, row.names = F){
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)}
    write.table(input, paste0(save_dir, "/", tsv_name), sep = "\t", row.names = row.names, col.names = T, quote = F)
}

#' Save ggplot Object as PDF
#'
#' This function saves a ggplot object as a PDF file with specified dimensions.
#' It tries to use Cairo PDF device if available, or falls back to standard PDF if not.
#'
#' @param input ggplot object to be saved
#' @param plot_name Name of the output PDF file
#' @param save_dir Directory where the PDF file will be saved
#' @param w Width of the PDF in inches
#' @param h Height of the PDF in inches
#' @return The file path (invisibly)
#' @examples
#' \dontrun{
#' # Save plot with default dimensions
#' save_plot(my_plot, "plot.pdf", "figures", w = 8, h = 6)
#' 
#' # Save plot with custom dimensions
#' save_plot(my_plot, "plot.pdf", "figures", w = 10, h = 8)
#' }
#' @export
save_plot <- function(input, plot_name, save_dir, w, h){
    # Create directory if it doesn't exist
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }
    
    # Construct file path
    file_path <- file.path(save_dir, plot_name)
    
    # Try using Cairo if available, otherwise use standard PDF
    tryCatch({
        if (requireNamespace("Cairo", quietly = TRUE)) {
            Cairo::CairoPDF(file = file_path, width = w, height = h)
        } else {
            pdf(file = file_path, width = w, height = h)
        }
        print(input)
        dev.off()
    }, error = function(e) {
        # If Cairo fails, try standard PDF
        message("Cairo PDF failed, using standard PDF device instead: ", e$message)
        pdf(file = file_path, width = w, height = h)
        print(input)
        dev.off()
    })
    
    # Return the file path invisibly
    invisible(file_path)
}

#' Save Expression Data from DESeq2 Object
#'
#' This function saves expression data from a DESeq2 object in various formats,
#' including the DESeq2 object itself, raw counts, normalized expression values,
#' and class labels for GSEA.
#'
#' @param dds DESeq2 object containing the expression data
#' @param save_dir Directory to save files. Default is the current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return Nothing, but saves the following files:
#'         - DESeq2 object as RDS
#'         - Raw counts as TSV
#'         - Variance-stabilized transformed data as TSV
#'         - CLS file with class labels
#' @examples
#' \dontrun{
#' # Save expression data with default settings
#' save_expression(dds)
#' 
#' # Save expression data with custom grouping and directory name
#' save_expression(dds, var = "Treatment", save_dir_name = "custom_results")
#' }
#' @export
save_expression <- function(dds, var, save_dir = getwd()){
    message("saving DESeq2 object & expressions...")
    counts <- assay(dds)
    data <- assay(vst(dds, blind = T))

    saveRDS(dds, file = paste0(save_dir, "/dds_qc.rds"))
    save_tsv(counts, tsv_name = "dds_counts.txt", save_dir = save_dir, row.names = T)
    save_tsv(data, tsv_name = "dds_vst.txt", save_dir = save_dir, row.names = T)
}