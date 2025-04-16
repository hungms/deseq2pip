#' Run Quality Control Pipeline
#'
#' This function runs the quality control portion of the RNA-seq analysis pipeline,
#' including filtering genes, generating QC plots, and saving expression data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param quantile Quantile threshold for filtering lowly expressed genes. Default is 0.05
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run QC pipeline with default settings
#' dds <- run_qc_pip(dds)
#' 
#' # Run QC pipeline with custom settings
#' dds <- run_qc_pip(dds, org = "mouse", remove_xy = TRUE,
#'                   save_dir_name = "custom_results")
#' }
#' @export
run_qc_pip <- function(
    dds,
    save_dir = getwd(), 
    org = "human",
    remove_xy = FALSE,
    remove_mt = FALSE,
    group_by = "Group1",
    quantile = 0.05,
    save_dir_name = "qc_results") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(org %in% c("mouse", "human"))
    stopifnot(is.logical(remove_xy))
    stopifnot(is.logical(remove_mt))
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(is.numeric(quantile) & quantile < 0.2)
    setwd(save_dir)
    
    # Remove XY genes if requested
    if(remove_xy){
        dds <- remove_xy_genes(dds, org = org)}
    
    # Remove MT genes if requested
    if(remove_mt){
        dds <- remove_mt_genes(dds, org = org)}
    
    # Filter lowly expressed genes
    dds <- remove_low_expression(dds, quantile = quantile, group_by = group_by, save_dir = save_dir, save_dir_name = save_dir_name)
    check_library(dds, save_dir = save_dir, save_dir_name = save_dir_name)

    # Save expression data
    save_expression(dds, group_by = group_by, save_dir = save_dir, save_dir_name = save_dir_name)
    
    return(dds)
}

#' Run Sample Distance Pipeline
#'
#' This function runs the distance analysis portion of the RNA-seq analysis pipeline,
#' including generating PCA and distance plots, and saving distance data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run distance pipeline with default settings
#' dds <- run_dist_pip(dds)
#' 
#' # Run distance pipeline with custom settings and batch correction
#' dds <- run_dist_pip(dds, batch = "Batch", save_dir_name = "qc_results")
#' }
#' @export
run_dist_pip <- function(
    dds,
    save_dir = getwd(), 
    group_by = "Group1",
    batch = NULL,
    pals = NULL,
    save_dir_name = "qc_results") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    if(length(batch) > 0){
        stopifnot(batch %in% colnames(colData(dds)))}
    if(length(pals) > 0){
        stopifnot(all(colData(dds)[[group_by]] %in% names(pals)))}
    setwd(save_dir)

    # batch correction
    message("Running batch correction...")
    vsd_default <- vst(dds, blind=FALSE)
    colnames(vsd_default) <- colnames(dds)
    
    # Initialize list with non-batch corrected data
    vsd_list <- list(default = vsd_default)
    if(length(batch) > 0){

        ## ComBat
        #dds_combat <- dds
        #assay(dds_combat) <- ComBat_seq(as.matrix(assay(dds)), batch=dds[[batch]])
        #vsd_combat <- vst(dds_combat, blind=FALSE)

        ## limma
        vsd_limma <- vsd_nobatch
        assay(vsd_limma) <- removeBatchEffect(as.matrix(assay(vsd_nobatch)), batch=vsd_nobatch[[batch]])
        vsd_list <- c(vsd_list, list(limma = vsd_limma))} #combat = vsd_combat, 

    # For each batch correction method
    message("Running PCA and distance analysis...")
    for (i in seq_along(vsd_list)) {
        select.vsd <- vsd_list[[i]]
        method_name <- names(vsd_list)[i]

        # Run PCA for additional nfcore groups if present
        nfcore.groups <- colnames(colData(select.vsd))[which(str_detect(colnames(colData(select.vsd)), "^Group"))]
        if (length(nfcore.groups) > 0) {
            for (j in seq_along(nfcore.groups)) {
                run_pca(
                    vsd = select.vsd, 
                    method = method_name, 
                    group_by = nfcore.groups[j], 
                    size = 4, 
                    save_dir = save_dir, 
                    save_dir_name = save_dir_name)
            }
        }
            
        # Run PCA and distance analysis
        run_pca(
            vsd = select.vsd, 
            method = method_name, 
            group_by = group_by, 
            pals = pals,
            size = 4, 
            save_dir = save_dir, 
            save_dir_name = save_dir_name)
        run_distance(
            vsd = select.vsd, 
            method = method_name, 
            save_data = T, 
            save_dir = save_dir, 
            save_dir_name = save_dir_name)
        
        # save vsd object
        save_tsv(assay(select.vsd), tsv_name = paste0(method_name, "_exprs.tsv"), save_dir = paste0(save_dir, "/", save_dir_name))
        }
    
    return(dds)
}

#' Run Differential Expression Pipeline
#'
#' This function runs the differential expression analysis portion of the RNA-seq pipeline,
#' including running DESeq2 and generating volcano plots.
#'
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return A data frame containing differential expression results for all comparisons
#' @examples
#' \dontrun{
#' # Run differential expression pipeline
#' res <- run_diffexp_pip(dds)
#' 
#' # Run with custom grouping
#' res <- run_diffexp_pip(dds, group_by = "Treatment")
#' }
#' @export
run_diffexp_pip <- function(
    dds,
    org = "human",
    order = "pxfc",
    group_by = "Group1",
    save_dir = getwd()) {
    
    # Input validation
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(as.character(design(dds)) != "~ 1")

    # Run differential expression analysis
    res <- run_diffexp_list(dds, org = org, group_by = group_by, order = order, save_dir = save_dir)

    # Generate MA plots
    plot <- plot_ma_list(res, save_plot = TRUE, save_dir = save_dir)
    
    # Generate volcano plots
    plot <- plot_volcano_list(res, save_plot = TRUE, save_dir = save_dir)
    return(res)
}

#' Run Gene Set Enrichment Analysis Pipeline
#'
#' This function runs the gene set enrichment analysis portion of the RNA-seq pipeline,
#' including running GSEA and generating enrichment plots.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return A list of GSEA results for all comparisons
#' @examples
#' \dontrun{
#' # Run GSEA pipeline
#' gsea_results <- run_gsea_pip(res)
#' 
#' # Run with custom parameters
#' gsea_results <- run_gsea_pip(res, org = "mouse")
#' }
#' @export
run_gsea_pip <- function(
    res = NULL,
    org = "human",
    save_dir = getwd()) {
    
    # Input validation
    if(length(res) != 0){
        stopifnot(is.data.frame(res))
        stopifnot(c("gene", "rank") %in% colnames(res))}
    
    # Run GSEA
    gsea <- run_gsea_list(res = res, org = org, save_dir = save_dir)
    
    # Generate GSEA plots
    gsea.df <- read_gsea_tsv_list(merge = TRUE, data_dir = save_dir)
    plot <- plot_gsea_barplot_list(gsea.df, save_dir = save_dir)

    # Format for EnrichmentMap
    format_enrichmentmap(org = org, collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME"), data_dir = save_dir)
    
    return(gsea)
}