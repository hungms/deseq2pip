#' Run batch correction by Limma
#'
#' This function performs batch correction on the variance-stabilized transformed data using the Limma package.
#' It can help identify potential batch effects or outliers.
#'
#' @param dds A DESeq2 object containing the normalized gene expression data
#' @param method Method to use for batch correction. Default is "limma".
#' @param batch Column name in colData(vsd) to use for batch correction.
#' @param save_dir Directory to save the results. Default is the current working directory
#' @return A DESeq2 object containing the batch-corrected data
#' @export
run_batch_correction <- function(dds, method = "limma", batch, save_data = TRUE, save_dir = getwd()){

    dds <- validate_dds(dds)
    method <- validate_method(method)
    uncorrected <- vst(dds, blind=FALSE)
    colnames(uncorrected) <- colnames(dds)
    vsd_limma <- uncorrected

    if(method == "limma"){
        assay(vsd_limma) <- removeBatchEffect(as.matrix(assay(uncorrected)), batch=uncorrected[[batch]])}
    else{
        stopifnot("No method listed, please use 'limma' or 'combat'")}

    if(save_data){
        save_tsv(assay(vsd_limma), tsv_name = paste0("dds_exprs_", method, "_corrected_for_", batch, ".tsv"), save_dir = save_dir, row.names = TRUE)}

    return(vsd_limma)
}


#' Run Principal Component Analysis
#'
#' This function performs PCA on the variance-stabilized transformed data and generates a PCA plot.
#' It can visualize sample relationships and identify potential batch effects or outliers.
#'
#' @param vsd A DESeq2 object containing the normalized gene expression data
#' @param group_by Column name in colData(vsd) to group by. 
#' @param shape Column name in colData(vsd) to use for shape in the PCA plot. Default is NULL
#' @param size Size of points in the PCA plot. Default is 4
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param save_prefix Prefix for the save file. Default is NULL
#' @param save_data Logical. If TRUE, saves PCA results to TSV. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the PCA plot to PDF. Default is TRUE
#' @param save_dir Directory to save the results. Default is the current working directory
#' @return A ggplot object showing the PCA plot
#' @examples
#' \dontrun{
#' # Basic PCA plot
#' p <- run_pca(vsd)
#' 
#' # PCA plot with shape and save results
#' p <- run_pca(vsd, shape = "Batch", save_data = TRUE)
#' 
#' # PCA plot with custom colors
#' p <- run_pca(vsd, pals = c("red", "blue", "green"))
#' }
#' @export
run_pca <- function(vsd, group_by, shape = NULL, pals = NULL, size = 4, save_prefix = NULL, save_data = TRUE, save_plot = TRUE, save_dir = getwd()){

    # validations
    dds <- validate_dds_group_by(vsd, group_by)
    pals <- validate_pals(vsd, group_by, pals)
    shape <- validate_shape(vsd, group_by, shape)
    validate_paths(save_dir)

    # add shape if provided
    if(!is.null(shape)){
        group_by <- c(group_by, shape)}
    
    # run pca
    pcadf <- plotPCA(vsd, intgroup=group_by, returnData=TRUE)
    
    percentVar <- round(100 * attr(pcadf, "percentVar"))
    p <- ggplot(pcadf, aes_string("PC1", "PC2", color=group_by[1])) +
        ggalt::geom_encircle(aes_string(fill = group_by[1]), alpha = 0.3) +
        geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
        geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
        ggprism::theme_prism(border = T) +
        labs(x = paste0("PC1 (", round(percentVar[1], 1), "% variance)"), y = paste0("PC2 (", round(percentVar[2], 1), "% variance)"))

    # add pals if provided
    if(!is.null(pals)) {
        p <- p + scale_color_manual(values = pals) + scale_fill_manual(values = pals)}

    # add shape if provided
    if(!is.null(shape)){
        p <- p + 
            geom_point(aes_string(shape = group_by[2]), size = size) +
            guides(fill = guide_legend(title = ""), color = guide_legend(title = ""), shape = guide_legend(title = ""))}
    else{
        p <- p + 
            geom_point(size = size) +
            guides(fill = guide_legend(title = ""), color = guide_legend(title = ""))}

    # save data
    if(save_data){
        save_tsv(pcadf, tsv_name = paste0("pca_", paste0(as.character(group_by, save_prefix), collapse = "_"), ".tsv"), save_dir = save_dir, row.names = TRUE)}
    
    # save plot
    if(save_plot){
        save_plot(p, plot_name = paste0("pca_", paste0(as.character(group_by, save_prefix), collapse = "_"), ".pdf"), save_dir = save_dir, w = 7, h = 5)}
    
    print(p)
    return(p)
}

#' Calculate and Plot Sample Distances
#'
#' This function calculates the Euclidean distance between samples based on their expression profiles
#' and generates a heatmap visualization. It helps identify sample relationships and potential outliers.
#'
#' @param vsd A DESeq2 object containing the normalized gene expression data
#' @param save_prefix Prefix for the save file. Default is NULL
#' @param save_data Logical. If TRUE, saves distance matrix to TSV. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the distance heatmap to PDF. Default is TRUE
#' @param save_dir Directory to save the results. Default is the current working directory
#' @param ... Additional arguments passed to pheatmap()
#' @return A pheatmap object showing the sample distances
#' @examples
#' \dontrun{
#' # Generate and display heatmap
#' p <- run_distance(vsd)
#' 
#' # Generate and save heatmap with custom parameters
#' p <- run_distance(vsd, save_plot = TRUE, save_dir_name = "custom_results")
#' }
#' @export
run_distance <- function(vsd, save_prefix = NULL, save_data = TRUE, save_plot = TRUE, save_dir = getwd(),  ...){

    # validations
    vsd <- validate_dds(vsd)
    validate_paths(save_dir)

    # format data
    dist <- dist(t(assay(vsd)))
    dist_mat <- as.matrix(dist)
    rownames(dist_mat) <- paste(colnames(vsd))
    colnames(dist_mat) <- paste(colnames(vsd))
    colors <- colorRampPalette(brewer.pal(9, "RdBu"))(255)

    # run pheatmap
    p <- pheatmap(dist_mat,
        show_column_names = FALSE,
        clustering_distance_rows=dist,
        clustering_distance_cols=dist,
        heatmap_legend_param = list(title = "Euclidean\nDistance"),
        col=colors,
        ...)

    if(save_data){
        save_tsv(dist_mat, tsv_name = paste0(paste0("euclidean_distance", as.character(save_prefix), collapse = "_"), ".tsv"), save_dir = save_dir, row.names = TRUE)}
    
    if(save_plot){
        wscale <- ncol(vsd)
        hscale <- max(nchar(colnames(vsd)))
        save_plot(p, plot_name = paste0(paste0("euclidean_distance", as.character(save_prefix), collapse = "_"), ".pdf"), save_dir = save_dir, w = wscale*0.2+hscale*0.3, h = wscale*0.2+hscale*0.2)}

    print(p)
    return(p)
}

#' Run Sample Distance Pipeline
#'
#' This function runs the distance analysis portion of the RNA-seq analysis pipeline,
#' including generating PCA and distance plots, and saving distance data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
#' @param save_dir Directory where all output files will be saved. Default is current working directory
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
    group_by,
    pals = NULL,
    batch = NULL,
    save_dir = getwd(), 
    save_dir_name = "qc_results") {
    
    # validations
    validate_paths(save_dir)
    dds <- validate_dds_group_by(dds, group_by)
    dds <- validate_batch(dds, batch)
    pals <- validate_pals(dds, group_by, pals)

    # create save directory
    qc_save_dir <- paste0(save_dir, "/", save_dir_name)
    dir.create(qc_save_dir, showWarnings = FALSE, recursive = TRUE)

    # batch correction
    message("Running batch correction...")
    vsd_default <- vst(dds, blind=FALSE)
    colnames(vsd_default) <- colnames(dds)
    
    # Initialize list with non-batch corrected data
    vsd_list <- list(default = vsd_default)
    if(!is.null(batch)){
        vsd_limma <- run_batch_correction(dds, method = "limma", batch, save_data = TRUE, save_dir = qc_save_dir)
        vsd_list <- c(vsd_list, list(limma = vsd_limma))}

    # For each batch correction method
    message("Running PCA and distance analysis...")
    for (i in seq_along(vsd_list)) {
        select.vsd <- vsd_list[[i]]
        method_name <- names(vsd_list)[i]
        if(method_name == "default"){save_prefix <- NULL}
        else{save_prefix <- paste0(method_name, "_corrected_for_", batch)}

        # Run PCA for additional nfcore groups if present
        nfcore.groups <- colnames(colData(select.vsd))[which(str_detect(colnames(colData(select.vsd)), "^Group"))]
        if (length(nfcore.groups) > 0) {
            for (j in seq_along(nfcore.groups)) {
                run_pca(
                    vsd = select.vsd, 
                    group_by = nfcore.groups[j], 
                    size = 4, 
                    save_prefix = save_prefix,
                    save_dir = qc_save_dir)
            }
        }
            
        # Run PCA and distance analysis
        run_pca(
            vsd = select.vsd, 
            group_by = group_by, 
            pals = pals,
            size = 4, 
            save_prefix = save_prefix,
            save_dir = qc_save_dir)
        run_distance(
            vsd = select.vsd, 
            save_prefix = save_prefix,
            save_data = T, 
            save_dir = qc_save_dir)
        }
    
    return(dds)
}


