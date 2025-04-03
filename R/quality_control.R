#' Remove XY Chromosome Genes from DESeq2 Object
#'
#' This function removes genes located on the X and Y chromosomes from a DESeq2 object.
#' It uses the Ensembl database to identify genes on these chromosomes.
#'
#' @param dds A DESeq2 object containing the gene expression data
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @return A filtered DESeq2 object with XY chromosome genes removed
#' @examples
#' \dontrun{
#' # For human data
#' dds_filtered <- remove_xy_genes(dds, org = "human")
#' 
#' # For mouse data
#' dds_filtered <- remove_xy_genes(dds, org = "mouse")
#' }
#' @export
remove_xy_genes <- function(dds, org = "human", ...){
    stopifnot(org %in% c("mouse", "human"))
    stopifnot("gene" %in% colnames(rowData(dds)))
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
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @return A filtered DESeq2 object with mitochondrial genes removed
#' @examples
#' \dontrun{
#' # For human data
#' dds_filtered <- remove_mt_genes(dds, org = "human")
#' 
#' # For mouse data
#' dds_filtered <- remove_mt_genes(dds, org = "mouse")
#' }
#' @export
remove_mt_genes <- function(dds, org = "human", ...){
    stopifnot(org %in% c("mouse", "human"))
    stopifnot("gene" %in% colnames(rowData(dds)))
    mt.genes <- get_mt_genes(org = org, ...)
    count <- count(rowData(dds)$gene %in% mt.genes)
    message(paste0("Removing ", count, " MT genes out of ", nrow(dds), " total genes..."))
    dds <- dds[which(!rowData(dds)$gene %in% mt.genes),]
    return(dds)
}

#' Remove Lowly Expressed Genes from DESeq2 Object
#'
#' This function removes genes with low expression levels from a DESeq2 object based on a quantile threshold.
#' It generates a density plot of expression levels and removes genes below the specified quantile threshold.
#'
#' @param dds A DESeq2 object containing the gene expression data
#' @param experiment Name of the experiment. Default is NULL
#' @param quantile The quantile threshold for low expression (0-1). Default is 0.05 (5th percentile)
#' @param group_by Column name in colData(dds) to group by. Default is "Group1"
#' @param save_plot Logical. If TRUE, saves the expression density plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is the current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return A filtered DESeq2 object with lowly expressed genes removed
#' @examples
#' \dontrun{
#' # Remove bottom 5% of genes
#' dds_filtered <- remove_low_expression(dds, quantile = 0.05)
#' 
#' # Remove bottom 10% of genes and save plot
#' dds_filtered <- remove_low_expression(dds, quantile = 0.1, save_plot = TRUE, save_dir_name = "custom_results")
#' }
#' @export
remove_low_expression <- function(dds, experiment = NULL, quantile = 0.05, group_by = "Group1", save_plot = TRUE, save_dir = getwd(), save_dir_name = "qc_results"){
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
        save_plot(p, experiment = experiment, plot_name = "low_expression.pdf", save_dir = paste0(save_dir, "/", save_dir_name, "/"), w=8, h=4)}

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
#' @param experiment Name of the experiment. Default is NULL
#' @param save_plot Logical. If TRUE, saves the library size plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is the current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
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
check_library <- function(dds, experiment = NULL, save_plot = TRUE, save_dir = getwd(), save_dir_name = "qc_results"){
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
        save_plot(p, experiment = experiment, plot_name = "library_size_distribution.pdf", save_dir = paste0(save_dir, "/", save_dir_name, "/"), w=0.4*wscale, h=3+0.2*hscale)}
    }

#' Run Principal Component Analysis
#'
#' This function performs PCA on the variance-stabilized transformed data and generates a PCA plot.
#' It can visualize sample relationships and identify potential batch effects or outliers.
#'
#' @param vsd A DESeq2 object containing the normalized gene expression data
#' @param experiment Name of the experiment. Default is NULL
#' @param group_by Column name in colData(vsd) to group by. Default is "Group1"
#' @param shape Column name in colData(vsd) to use for shape in the PCA plot. Default is NULL
#' @param size Size of points in the PCA plot. Default is 4
#' @param cols Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param save_data Logical. If TRUE, saves PCA results to TSV. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the PCA plot to PDF. Default is TRUE
#' @param save_dir Directory to save the results. Default is the current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return A ggplot object showing the PCA plot
#' @examples
#' \dontrun{
#' # Basic PCA plot
#' p <- run_pca(vsd)
#' 
#' # PCA plot with shape and save results
#' p <- run_pca(vsd, shape = "Batch", save_data = TRUE, save_dir_name = "custom_results")
#' 
#' # PCA plot with custom colors
#' p <- run_pca(vsd, cols = c("red", "blue", "green"))
#' }
#' @export
run_pca <- function(vsd, experiment = NULL, group_by = "Group1", shape = NULL, size = 4, cols = NULL, save_data = TRUE, save_plot = TRUE, save_dir = getwd(), save_dir_name = "qc_results"){
    stopifnot(length(group_by) == 1)

    if(length(shape) > 0){
        group_by <- c(group_by, shape)}
    
    pcadf <- plotPCA(vsd, intgroup=group_by, returnData=TRUE)
    p <- ggplot(pcadf, aes_string("PC1", "PC2", color=group_by[1])) +
        ggalt::geom_encircle(aes_string(fill = group_by[1]), alpha = 0.3) +
        umap_aes() +
        theme_text()

    if(!is.null(cols)) {
        p <- p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
    }

    if(length(shape) > 0){
        p <- p + 
            geom_point(aes_string(shape = group_by[2]), size = size) +
            guides(fill = guide_legend(title = ""), color = guide_legend(title = ""), shape = guide_legend(title = ""))}
    else{
        p <- p + 
            geom_point(size = size) +
            guides(fill = guide_legend(title = ""), color = guide_legend(title = ""))}
    
    print(p)

    if(save_data){
        save_tsv(pcadf, experiment = experiment, tsv_name = paste0("pca_", paste0(group_by, collapse = "_"), ".tsv"), save_dir = paste0(save_dir, "/", save_dir_name, "/"), row.names = TRUE)}
    if(save_plot){
        save_plot(p, experiment = experiment, plot_name = paste0("pca_", paste0(group_by, collapse = "_"), ".pdf"), save_dir = paste0(save_dir, "/", save_dir_name, "/"), w = 7, h = 5)}
    
    return(p)
}

#' Calculate and Plot Sample Distances
#'
#' This function calculates the Euclidean distance between samples based on their expression profiles
#' and generates a heatmap visualization. It helps identify sample relationships and potential outliers.
#'
#' @param vsd A DESeq2 object containing the normalized gene expression data
#' @param experiment Name of the experiment. Default is NULL
#' @param save_data Logical. If TRUE, saves distance matrix to TSV. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the distance heatmap to PDF. Default is TRUE
#' @param save_dir Directory to save the results. Default is the current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
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
run_distance <- function(vsd, experiment = NULL, save_data = TRUE, save_plot = TRUE, save_dir = getwd(), save_dir_name = "qc_results", ...){
    dist <- dist(t(assay(vsd)))
    dist_mat <- as.matrix(dist)
    rownames(dist_mat) <- paste(colnames(vsd))
    colnames(dist_mat) <- paste(colnames(vsd))
    colors <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
    p <- pheatmap(dist_mat,
            show_column_names = FALSE,
            clustering_distance_rows=dist,
            clustering_distance_cols=dist,
            heatmap_legend_param = list(title = "Euclidean\nDistance"),
            col=colors,
            ...)
    print(p)

    if(save_data){
        save_tsv(dist_mat, experiment = experiment, tsv_name = paste0("euclidean_distance.tsv"), save_dir = paste0(save_dir, "/", save_dir_name, "/"), row.names = TRUE)}
    
    if(save_plot){
        wscale <- ncol(vsd)
        hscale <- max(nchar(colnames(vsd)))
        save_plot(p, experiment = experiment, plot_name = paste0("euclidean_distance_heatmap.pdf"), save_dir =  paste0(save_dir, "/", save_dir_name, "/"), w = wscale*0.2+hscale*0.3, h = wscale*0.2+hscale*0.2)}

    return(p)
}


