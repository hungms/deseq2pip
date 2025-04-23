#' Remove XY Chromosome Genes from DESeq2 Object
#'
#' This function removes genes located on the X and Y chromosomes from a DESeq2 object.
#' It uses the Ensembl database to identify genes on these chromosomes.
#'
#' @param dds A DESeq2 object containing the gene expression data
#' @param org The organism to use, either "human" or "mouse".
#' @return A filtered DESeq2 object with XY chromosome genes removed
#' @examples
#' \dontrun{
#' # For human data
#' dds_filtered <- remove_xy_genes(dds, org)
#' 
#' # For mouse data
#' dds_filtered <- remove_xy_genes(dds, org = "mouse")
#' }
#' @export
remove_xy_genes <- function(dds, org, ...){
    dds <- validate_dds(dds)
    org <- validate_org(org)
    xy.genes <- get_xy_genes(org = org, ...)
    count <- count(rowData(dds)$gene %in% xy.genes)
    message(paste0("Removing ", count, " XY genes out of ", nrow(dds), " total genes..."))
    dds <- dds[which(!rowData(dds)$gene %in% xy.genes),]
    return(dds)
}

#' Remove Mitochondrial Genes from DESeq2 Object
#'
#' This function removes genes encoded in the mitochondrial genome from a DESeq2 object.
#' It uses the Ensembl database to identify mitochondrial genes.
#'
#' @param dds A DESeq2 object containing the gene expression data
#' @param org The organism to use, either "human" or "mouse".
#' @return A filtered DESeq2 object with mitochondrial genes removed
#' @examples
#' \dontrun{
#' # For human data
#' dds_filtered <- remove_mt_genes(dds, org)
#' 
#' # For mouse data
#' dds_filtered <- remove_mt_genes(dds, org = "mouse")
#' }
#' @export
remove_mt_genes <- function(dds, org, ...){
    dds <- validate_dds(dds)
    org <- validate_org(org)
    mt.genes <- get_mt_genes(org = org, ...)
    count <- count(rowData(dds)$gene %in% mt.genes)
    message(paste0("Removing ", count, " MT genes out of ", nrow(dds), " total genes..."))
    dds <- dds[which(!rowData(dds)$gene %in% mt.genes),]
    return(dds)
}

#' Remove Genes with Low Expression
#'
#' This function filters out genes with expression values below a specified quantile threshold
#' in at least a minimum number of replicates per condition.
#'
#' @param dds A DESeq2 object containing the gene expression data
#' @param quantile Quantile threshold for filtering. Default is 0.05 (lowest 5%)
#' @param group_by Column name in colData(dds) to use for defining conditions. Default is "Group1"
#' @param save_plot Logical. If TRUE, saves the expression distribution plot. Default is TRUE
#' @param save_dir Directory to save the plot. Default is the current working directory
#' @return A filtered DESeq2 object with low-expression genes removed
#' @examples
#' \dontrun{
#' # Remove bottom 5% of genes
#' dds_filtered <- remove_low_expression(dds, quantile = 0.05)
#' 
#' # Remove bottom 10% of genes and save plot
#' dds_filtered <- remove_low_expression(dds, quantile = 0.1, save_plot = TRUE)
#' }
#' @export
remove_low_expression <- function(dds, group_by, quantile = 0.05, save_plot = TRUE, save_dir = getwd()){

    # validations
    dds <- validate_dds_group_by(dds, group_by)
    quantile <- validate_quantile(quantile)
    validate_paths(save_dir)
    
    # run filter
    message("Filtering genes with low expressions...")
    vsd <- vst(dds, blind=FALSE)
    colnames(vsd) <- colnames(dds)
    
    vsd.mean <- rowMeans(assay(vsd)) %>%
        as.numeric(.) %>%
        as.data.frame(.)
    threshold <- quantile(vsd.mean[[1]], quantile)

    p <- vsd.mean %>%
            ggplot(aes(x = .)) +
            geom_density(size = 0.6, color = "grey40") +
            geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
            xlab("Expression") +
            ylab("Density") +
            theme_line() +
            theme_text()
    print(p)

    if(save_plot){
        save_plot(p, plot_name = "low_expression.pdf", save_dir = save_dir, w=8, h=4)}

    min_rep <- min(table(colData(dds)[[group_by]]))
    keep <- rowSums(assay(vsd) >= threshold) >= min_rep
    dds <- dds[keep,]
    return(dds)}

#' Check Library Size Distribution
#'
#' This function generates a boxplot showing the distribution of expression values across samples.
#' It helps identify potential outliers or quality issues in the sequencing libraries.
#'
#' @param dds A DESeq2 object containing the gene expression data
#' @param save_plot Logical. If TRUE, saves the library size plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is the current working directory
#' @return A ggplot object showing the library size distribution
#' @examples
#' \dontrun{
#' # Generate and display plot
#' p <- check_library(dds)
#' 
#' # Generate and save plot
#' p <- check_library(dds, save_plot = TRUE, save_dir_name = "custom_results")
#' }
#' @export
check_library <- function(dds, save_plot = TRUE, save_dir = getwd()){
    vsd <- vst(dds, blind=FALSE) 
    assay <- assay(vsd) %>% as.data.frame(.)
    colnames(vsd) <- colnames(dds)
    p <- assay %>%
        pivot_longer(everything(), names_to = "samples", values_to = "exprs") %>%
        ggplot(aes(x = samples, y = exprs)) +
        geom_boxplot(fill = "grey", width = 0.75) +
        theme_border() +
        theme_text() +
        xlab(NULL) +
        ylab("Expression") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    print(p)

    if(save_plot){
        wscale <- ncol(dds)
        hscale <- max(nchar(colnames(dds)))
        save_plot(p, plot_name = "library_size_distribution.pdf", save_dir = save_dir, w=0.4*wscale, h=3+0.2*hscale)}
    }

#' Run Quality Control Pipeline
#'
#' This function runs the quality control portion of the RNA-seq analysis pipeline,
#' including filtering genes, generating QC plots, and saving expression data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse".
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is TRUE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is TRUE
#' @param group_by Column name in colData(dds) to use for grouping.
#' @param quantile Quantile threshold for filtering lowly expressed genes. Default is 0.05
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run QC pipeline with default settings
#' dds <- run_qc_pip(dds)
#' 
#' # Run QC pipeline with custom settings
#' dds <- run_qc_pip(dds, org = "mouse", remove_xy = TRUE, save_dir_name = "custom_results")
#' }
#' @export
run_qc_pip <- function(
    dds,
    org,
    remove_xy = TRUE,
    remove_mt = TRUE,
    group_by,
    quantile = 0.05,
    save_dir = getwd(), 
    save_dir_name = "qc_results") {
    
    # validations
    dds <- validate_dds_group_by(dds, group_by)
    org <- validate_org(org)
    quantile <- validate_quantile(quantile)
    validate_logical(c(remove_xy, remove_mt))
    validate_paths(save_dir)

    # create save directory
    qc_save_dir <- paste0(save_dir, "/", save_dir_name)
    dir.create(qc_save_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Remove XY genes if requested
    if(remove_xy){
        dds <- remove_xy_genes(dds, org = org)}

    # Remove MT genes if requested
    if(remove_mt){
        dds <- remove_mt_genes(dds, org = org)}
    
    # Filter lowly expressed genes
    dds <- remove_low_expression(dds, quantile = quantile, group_by = group_by, save_dir = qc_save_dir)
    check_library(dds, save_dir = qc_save_dir)

    # Save expression data
    save_expression(dds, group_by = group_by, save_dir = qc_save_dir)
    
    return(dds)
}