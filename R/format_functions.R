#' Summarize Gene Expression in DESeq2 Object
#'
#' This function aggregates expression values from multiple isoforms/peaks of the same gene
#' in a DESeq2 object into a single gene-level expression value. It creates a new
#' DESeq2 object with gene-level expression data.
#'
#' @param dds DESeq2 object containing isoform/peak-level expression data
#' @param gene_sym_col Column name in rowData(dds) containing gene symbols. Default is "Gene.Name"
#' @param ... Additional arguments passed to summarize_genes()
#' @return A new DESeq2 object containing gene-level expression data
#' @examples
#' \dontrun{
#' # Sum gene expression
#' gene_dds <- summarize_genes_dds(dds)
#' 
#' # Average gene expression
#' gene_dds <- summarize_genes_dds(dds, normalized = TRUE)
#' }
#' @export
summarize_genes_dds <- function(dds, gene_sym_col = "gene", ...){
    df <- assay(dds) %>% as.data.frame(.)
    gene_sym_vec <- rowData(dds)[[gene_sym_col]]
    df <- strpip::summarize_genes(df = df, gene_sym_vec = gene_sym_vec, ...)

    metadata <- as.data.frame(colData(dds))[colnames(df),]
    dds.new <- DESeqDataSetFromMatrix(
        countData = df,
        colData = metadata,
        design = design(dds))
    return(dds.new)
}

#' Import nfcore/rnaseq DESeq2 Object
#'
#' This function imports a DESeq2 object from nfcore RNA-seq pipeline output and adds gene names.
#'
#' @param rdata Path to the RData file containing the DESeq2 object
#' @param tx2gene Path to the tx2gene mapping file
#' @return A DESeq2 object with gene names added to rowData
#' @examples
#' \dontrun{
#' # Import DESeq2 object and add gene names
#' dds <- import_nfcore_rna("path/to/dds.RData", "path/to/tx2gene.tsv")
#' }
#' @export
import_nfcore_rna <- function(rdata, tx2gene){
    # Load DESeq2 object
    load(rdata)
    dds <- get("dds")

    # Read gene symbol mapping and ensure it's properly formatted
    genes <- read.table(tx2gene, sep = "\t", header = F) %>%
        dplyr::mutate(V3 = as.character(V3)) %>%  # Explicitly convert gene symbols to character
        dplyr::distinct(V2, V3) %>%
        tibble::column_to_rownames("V2")
    
    # Ensure rownames of dds are character vectors
    rownames(dds) <- as.character(rownames(dds))
    rowData(dds)$gene <-  as.character(genes[rownames(dds),])

    # Summarize genes
    dds.summarized <- summarize_genes_dds(dds, gene_sym_col = "gene", normalize = F)
    return(dds.summarized)
}

#' Import nfcore ATAC-seq DESeq2 Object
#'
#' This function imports a DESeq2 object from nfcore ATAC-seq pipeline output and adds gene names.
#'
#' @param rdata Path to the RData file containing the DESeq2 object
#' @param tx2gene Path to the tx2gene mapping file
#' @param dist.to.TSS Distance to TSS to consider as TSS. Default is 2000
#' @return A DESeq2 object with gene names added to rowData
#' @examples
#' \dontrun{
#' # Import DESeq2 object and add gene names
#' dds <- import_nfcore_atac("path/to/dds.RData", "path/to/tx2gene.tsv")
#' }
#' @export
import_nfcore_atac <- function(rdata, annotatePeaks, dist.to.TSS = 2000){

    # Load DESeq2 object
    load(rdata)
    dds <- get("dds")
    
    annotations <- read.table(annotatePeaks, sep = "\t", header = T, row.names = 1)
    annotations <- annotations[rownames(assay(dds)),]
    rowData(dds) <- cbind(rowData(dds), annotations)
    rowData(dds)$gene <- rowData(dds)$Gene.Name
    rowData(dds)$Gene.Name <- NULL
    rowData(dds)$TSS <- ifelse(rowData(dds)$gene != "" & rowData(dds)$Distance.to.TSS > -dist.to.TSS & rowData(dds)$Distance.to.TSS < dist.to.TSS, T, F)

    return(dds)
}

#' Get TSS peaks from DESeq2 object
#'
#' This function extracts TSS peaks from a DESeq2 object and removes duplicate genes.
#'
#' @param dds A DESeq2 object
#' @return A DESeq2 object containing only TSS peaks
#' @examples
#' \dontrun{
#' # Get TSS peaks
#' dds.tss <- getTSS(dds)
#' }
#' @export
getTSS <- function(dds){
    stopifnot("TSS" %in% colnames(rowData(dds)))
    dds.tss <- dds[which(rowData(dds)$TSS == TRUE), ]
    duplicated.genes <- rowData(dds.tss)$gene[duplicated(rowData(dds.tss)$gene)]
    dds.tss <- dds.tss[-c(which(rowData(dds.tss)$gene %in% duplicated.genes)),]
    rownames(dds.tss) <- rowData(dds.tss)$gene
    return(dds.tss)}

#' Import MSigDB Gene Sets
#'
#' This function imports pre-defined MSigDB gene sets for human or mouse organisms.
#'
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @return A data frame containing MSigDB gene sets
#' @examples
#' \dontrun{
#' # Import human MSigDB gene sets
#' msigdb_human <- import_msigdb("human")
#' 
#' # Import mouse MSigDB gene sets
#' msigdb_mouse <- import_msigdb("mouse")
#' }
#' @export
import_msigdb <- function(org) {
    files <- list.files(system.file("extdata", package = "deseq2pip"), full.names = T)
    files <- files[which(str_detect(files, paste0(org, "_msigdbr.tsv")))]
    f <- files[length(files)]
    msigdb <- read.table(gzfile(f), header = TRUE, sep = "\t")
    return(msigdb)}

#' Save Data Frame as TSV File
#'
#' This function saves a data frame as a tab-separated values (TSV) file.
#'
#' @param input Data frame to be saved
#' @param experiment Name of the experiment. Default is NULL
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
save_tsv <- function(input, experiment = NULL, tsv_name, save_dir, row.names = F){
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)}
    if(length(experiment) == 0){
        experiment <- ""}
    else{
        experiment <- paste0(experiment, "_")}
    write.table(input, paste0(save_dir, "/", experiment, tsv_name), sep = "\t", row.names = row.names, col.names = T, quote = F)}

#' Save ggplot Object as PDF
#'
#' This function saves a ggplot object as a PDF file with specified dimensions.
#' It tries to use Cairo PDF device if available, or falls back to standard PDF if not.
#'
#' @param input ggplot object to be saved
#' @param experiment Name of the experiment. Default is NULL
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
save_plot <- function(input, experiment = NULL, plot_name, save_dir, w, h){
    # Create directory if it doesn't exist
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }
    
    # Format experiment name
    if(length(experiment) == 0){
        experiment <- ""
    } else {
        experiment <- paste0(experiment, "_")
    }
    
    # Construct file path
    file_path <- file.path(save_dir, paste0(experiment, plot_name))
    
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
#' @param experiment Name of the experiment
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
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
#' save_expression(dds, experiment = "my_experiment")
#' 
#' # Save expression data with custom grouping and directory name
#' save_expression(dds, group_by = "Treatment", experiment = "my_experiment", save_dir_name = "custom_results")
#' }
#' @export
save_expression <- function(dds, experiment, group_by = "Group1", save_dir = getwd(), save_dir_name = "qc_results"){
    message("Saving DESeq2 object & expressions...")
    counts <- assay(dds)
    data <- assay(vst(dds, blind = T))

    group.lv <- rev(levels(dds[[group_by]]))
    ngroup <- length(group.lv)
    add_pad <- function(row, length) {
        c(row, rep("", length - length(row)))}
    row1 <- add_pad(c(ncol(dds), ngroup, 1), ncol(dds))
    row2 <- add_pad(group.lv, ncol(dds))
    row3 <- colnames(dds)

    class <- matrix("", nrow = 3, ncol = ncol(dds))
    class[1,] <- row1
    class[2,] <- row2
    class[3,] <- row3

    saveRDS(dds, file = paste0(save_dir, "/", save_dir_name, "/", experiment, "_dds_qc.rds"))
    save_tsv(counts, experiment = experiment, tsv_name = "expr.txt", save_dir = paste0(save_dir, "/", save_dir_name, "/"), row.names = T)
    save_tsv(data, experiment = experiment, tsv_name = "vst.txt", save_dir = paste0(save_dir, "/", save_dir_name, "/"), row.names = T)
    write.table(class, file = paste0(save_dir, "/", save_dir_name, "/", experiment, "_class.cls"), row.names = F, col.names = F, quote = F)
}