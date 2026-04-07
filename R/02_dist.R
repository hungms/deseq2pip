#' Run Sample Distance Pipeline

#' This function runs the distance analysis portion of the RNA-seq analysis pipeline,
#' including generating PCA and distance plots, and saving distance data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param var Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param design Design formula for the DESeq2 object.
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
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
    var,
    batch = NULL,
    pals = NULL,
    save_dir = getwd(), 
    save_dir_name = "qc_results") {
    
    # validations
    dds <- validate_var(dds, var)
    if(!is.null(batch)) dds <- validate_batch(dds, batch)
    if(!is.null(pals)) pals <- validate_pals(dds, var, pals)

    message("\n########################################################\nRunning Distance Analysis Pipeline\n########################################################")

    # create save directory
    qc_save_dir <- paste0(save_dir, "/", save_dir_name)
    dir.create(qc_save_dir, showWarnings = FALSE, recursive = TRUE)

    # batch correction
    message("performing variance stabilizing transformation (VST)...")
    uncorrected_vsd <- vst(dds, blind=FALSE)
    colnames(uncorrected_vsd) <- colnames(dds)
    
    # initialize list to store uncorrected and corrected expressions
    exprs_list <- list(uncorrected = uncorrected_vsd)

    # perform batch correction
    if(!is.null(batch)){
        corrected_list <- run_batch_correction(dds, batch, method = "limma", save_data = TRUE, save_dir = qc_save_dir)
        exprs_list <- c(exprs_list, corrected_list)}

    # for each uncorrected and corrected expressions
    for (i in seq_along(exprs_list)) {

        # setup expression object
        select_vsd <- exprs_list[[i]]
        save_prefix <- names(exprs_list)[i]
            
        # run PCA for main group, batch, and other nfcore groups if present
        nfcore_groups <- colnames(colData(select_vsd))[which(str_detect(colnames(colData(select_vsd)), "^Group"))]
        pca_var <- unique(c(var, batch, nfcore_groups))
        pca_var <- pca_var[which(pca_var %in% colnames(colData(select_vsd)))]
        run_pca(
            vsd = select_vsd, 
            var = pca_var, 
            pals = pals,
            size = 4, 
            save_prefix = save_prefix,
            save_dir = qc_save_dir)

        # run distance analysis
        run_distance(
            vsd = select_vsd, 
            save_prefix = save_prefix,
            save_data = T, 
            save_dir = qc_save_dir)
        }
    
    return(dds)
}


#' Run batch correction by Limma & Combat
#'
#' This function performs batch correction on the variance-stabilized transformed data using the Limma package.
#' It can help identify potential batch effects or outliers.
#'
#' @param dds A DESeq2 object containing the normalized gene expression data to be batch corrected
#' @param batch Column name in colData(vsd) to use for batch correction.
#' @param method Method to use for batch correction. Options are "limma" and "combat".
#' @param save_data Logical. If TRUE, saves the batch-corrected data to a TSV file. Default is TRUE
#' @param save_dir Directory to save the results. Default is the current working directory
#' @return A DESeq2 object containing the batch-corrected data
#' @export
run_batch_correction <- function(dds, batch, method = c("limma", "combat"), save_data = TRUE, save_dir = getwd()) {

    # validations
    dds <- validate_dds(dds)
    method <- validate_method(method)
    dds <- validate_batch(dds, batch)
    validate_logical(save_data)

    uncorrected <- vst(dds, blind=FALSE)
    colnames(uncorrected) <- colnames(dds)

    # run limma batch correction
    corrected_vsd <- list()
    if("limma" %in% method){
        message("performing limma batch correction...")
        limma_vsd <- uncorrected
        assay(limma_vsd) <- removeBatchEffect(as.matrix(assay(limma_vsd)), batch=limma_vsd[[batch]])
        corrected_vsd[["limma"]] <- limma_vsd}

    # run combat batch correction
    if("combat" %in% method){
        message("performing combat batch correction...")
        combat_vsd <- uncorrected
        assay(combat_vsd) <- sva::ComBat(as.matrix(assay(combat_vsd)), batch=combat_vsd[[batch]])
        corrected_vsd[["combat"]] <- combat_vsd}

    # save data
    if(save_data){
        for(i in seq_along(corrected_vsd)){
            save_tsv(assay(corrected_vsd[[i]]), tsv_name = paste0("dds_exprs_", names(corrected_vsd)[i], ".tsv"), save_dir = save_dir, row.names = TRUE)}
    }

    return(corrected_vsd)
}


#' Run Principal Component Analysis
#'
#' This function performs PCA on the variance-stabilized transformed data and generates a PCA plot.
#' It can visualize sample relationships and identify potential batch effects or outliers.
#'
#' @param vsd A DESeq2 object containing the normalized gene expression data
#' @param var Column names in colData(vsd) to group by. 
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
run_pca <- function(vsd, var, pals = NULL, size = 4, save_prefix = NULL, save_data = TRUE, save_plot = TRUE, save_dir = getwd()){

    # validations
    vsd <- validate_var(vsd, var)
    if(!is.null(pals) && length(var) == 1) pals <- validate_pals(vsd, var, pals)
    validate_logical(c(save_data, save_plot))

    # run pca
    message("running PCA analysis...")
    pcadf <- plotPCA(vsd, intgroup=var, returnData=TRUE)

    # plot pca
    pca_list <- list()
    for(i in seq_along(var)){

        message("\t- generating PCA plot for ", var[i], "...")
        percentVar <- round(100 * attr(pcadf, "percentVar"))
        p <- ggplot(pcadf, aes_string("PC1", "PC2", color=var[i])) +
            #ggalt::geom_encircle(aes_string(fill = var[i]), alpha = 0.3) +
            geom_point(size = size) +
            geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
            geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
            ggprism::theme_prism(border = T) +
            labs(x = paste0("PC1 (", round(percentVar[1], 1), "% variance)"), y = paste0("PC2 (", round(percentVar[2], 1), "% variance)"))

        # add pals if provided
        if(!is.null(pals)) {
            if(all(names(pals) %in% pcadf[[var[i]]])) {
                p <- p + scale_color_manual(values = pals) + scale_fill_manual(values = pals)}
            else{
                message("names(pals) do not match the levels of ", var[i], "...")}}

        # save plot
        if(save_plot){
            save_plot(p, plot_name = paste0("pca_", paste0(as.character(var[i]), "_", save_prefix), ".pdf"), save_dir = save_dir, w = 7, h = 5)}
        
        print(p)
        pca_list[[i]] <- p}

    # save data
    if(save_data){
        message("saving PCA results to TSV...")
        save_tsv(pcadf, tsv_name = paste0("pca_", save_prefix, ".tsv"), save_dir = save_dir, row.names = TRUE)}
    return(pca_list)
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
    validate_logical(c(save_data, save_plot))

    # calculate euclidean distance
    message("calculating euclidean distance between samples...")
    dist <- dist(t(assay(vsd)))
    dist_mat <- as.matrix(dist)
    rownames(dist_mat) <- paste(colnames(vsd))
    colnames(dist_mat) <- paste(colnames(vsd))
    colors <- colorRampPalette(brewer.pal(9, "RdBu"))(255)

    # plot distance heatmap
    message("generating sample distance heatmap...")
    p <- pheatmap(dist_mat,
        show_column_names = FALSE,
        clustering_distance_rows=dist,
        clustering_distance_cols=dist,
        heatmap_legend_param = list(title = "Euclidean\nDistance"),
        col=colors,
        ...)

    # save data
    if(save_data){
        message("saving distance matrix to TSV...")
        save_tsv(dist_mat, tsv_name = paste0(paste0("euclidean_distance_", as.character(save_prefix), collapse = "_"), ".tsv"), save_dir = save_dir, row.names = TRUE)}
    
    # save plot
    if(save_plot){
        wscale <- ncol(vsd)
        hscale <- max(nchar(colnames(vsd)))
        save_plot(p, plot_name = paste0(paste0("euclidean_distance_", as.character(save_prefix), collapse = "_"), ".pdf"), save_dir = save_dir, w = wscale*0.2+hscale*0.3, h = wscale*0.2+hscale*0.2)}

    print(p)
    return(p)
}

