#' Run Quality Control Pipeline
#'
#' This function runs the quality control portion of the RNA-seq analysis pipeline,
#' including filtering genes, generating QC plots, and saving expression data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param experiment Name of the experiment
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param quantile Quantile threshold for filtering lowly expressed genes. Default is 0.05
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' # Run QC pipeline with default settings
#' dds <- run_qc_pip(dds, experiment = "my_experiment")
#' 
#' # Run QC pipeline with custom settings
#' dds <- run_qc_pip(dds, experiment = "my_experiment",
#'                   org = "mouse", remove_xy = TRUE,
#'                   save_dir_name = "custom_results")
#' @export
run_qc_pip <- function(
    dds,
    experiment,
    save_dir = getwd(), 
    org = "human",
    remove_xy = F,
    remove_mt = F,
    group_by = "Group1",
    quantile = 0.05,
    save_dir_name = "qc_results") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(length(experiment) == 1 & is.character(experiment))
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

    # Run PCA and distance analysis
    run_distance(dds, save_data = T, save_dir = save_dir, save_dir_name = save_dir_name)
    run_pca(dds, group_by = group_by, size = 4, save_dir = save_dir, save_dir_name = save_dir_name)
    
    # Run PCA for additional nfcore groups if present
    nfcore.groups <- colnames(colData(dds))[which(str_detect(colnames(colData(dds)), "^Group"))]
    if(length(nfcore.groups) > 0){
        for(i in seq_along(nfcore.groups)){
            run_pca(dds, group_by = nfcore.groups[i], size = 4, save_dir = save_dir, save_dir_name = save_dir_name)}}

    # Save expression data
    save_expression(dds, group_by = group_by, experiment = experiment, save_dir = save_dir, save_dir_name = save_dir_name)
    
    return(dds)
}

#' Run Differential Expression Pipeline
#'
#' This function runs the differential expression analysis portion of the RNA-seq pipeline,
#' including running DESeq2 and generating volcano plots.
#'
#' @param dds DESeq2 object containing the expression data
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return A data frame containing differential expression results for all comparisons
#' @examples
#' # Run differential expression pipeline
#' res <- run_diffexp_pip(dds)
#' 
#' # Run with custom grouping
#' res <- run_diffexp_pip(dds, group_by = "Treatment")
#' @export
run_diffexp_pip <- function(
    dds,
    org = "human",
    group_by = "Group1",
    save_dir = getwd()) {
    
    # Input validation
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(as.character(design(dds)) != "~ 1")
    
    # Run differential expression analysis
    res <- run_diffexp_list(dds, org = org, group_by = group_by, save_dir = save_dir)
    
    # Generate volcano plots
    plot <- plot_volcano_list(res, save_plot = T, save_dir = save_dir)
    
    return(res)
}

#' Run Gene Set Enrichment Analysis Pipeline
#'
#' This function runs the gene set enrichment analysis portion of the RNA-seq pipeline,
#' including running GSEA and generating enrichment plots.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param order Column name to use for ranking genes. Default is "rank"
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return A list of GSEA results for all comparisons
#' @examples
#' # Run GSEA pipeline
#' gsea_results <- run_gsea_pip(res)
#' 
#' # Run with custom parameters
#' gsea_results <- run_gsea_pip(res, org = "mouse")
#' @export
run_gsea_pip <- function(
    res = NULL,
    org = "human",
    order = "rank",
    save_dir = getwd()) {
    
    # Input validation
    if(length(res) != 0){
        stopifnot(is.data.frame(res))
        stopifnot(c("padj", "gene", "log2FoldChange", "comparison") %in% colnames(res))
        stopifnot(order %in% c("rank", "log2FoldChange", "pvalue", "padj"))}
    
    # Run GSEA
    gsea <- run_gsea_list(res = res, org = org, order = order, save_dir = save_dir)
    
    # Generate GSEA plots
    gsea.df <- read_gsea_tsv_list(merge = T, data_dir = save_dir)
    plot <- plot_gsea_barplot_list(gsea.df, save_dir = save_dir)

    # Format for EnrichmentMap
    format_enrichmentmap(collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME"), data_dir = save_dir)
    
    return(gsea)
}

#' Run Complete DESeq2 Pipeline
#'
#' This function runs a complete RNA-seq analysis pipeline including quality control,
#' differential expression analysis, and gene set enrichment analysis for all possible
#' comparisons in the experiment. It generates various plots and saves results in
#' organized directories.
#'
#' @param dds DESeq2 object containing the expression data
#' @param experiment Name of the experiment
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param quantile Quantile threshold for filtering lowly expressed genes. Default is 0.05
#' @param order Column name to use for ranking genes. Default is "rank"
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' # Run complete pipeline with default settings
#' dds <- run_deseq2_pip(dds, experiment = "my_experiment")
#' 
#' # Run pipeline with custom settings
#' dds <- run_deseq2_pip(dds, experiment = "my_experiment",
#'                     org = "mouse", remove_xy = TRUE,
#'                     group_by = "Treatment",
#'                     save_dir_name = "custom_results")
#' @export
run_deseq2_pip <- function(
    dds,
    experiment,
    save_dir = getwd(), 
    org = "human",
    remove_xy = F,
    remove_mt = F,
    group_by = "Group1",
    quantile = 0.05,
    order = "rank",
    save_dir_name = "qc_results") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(length(experiment) == 1 & is.character(experiment))
    stopifnot(org %in% c("mouse", "human"))
    stopifnot(is.logical(remove_xy))
    stopifnot(is.logical(remove_mt))
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(is.numeric(quantile) & quantile < 0.2)
    stopifnot(order %in% c("rank", "log2FoldChange", "pvalue", "padj"))
    stopifnot(as.character(design(dds)) != "~ 1")
    
    # Run quality control pipeline
    dds <- run_qc_pip(
        dds = dds,
        experiment = experiment,
        save_dir = save_dir,
        org = org,
        remove_xy = remove_xy,
        remove_mt = remove_mt,
        group_by = group_by,
        quantile = quantile,
        save_dir_name = save_dir_name
    )
    
    # Run differential expression pipeline
    res <- run_diffexp_pip(
        dds = dds,
        org = org,
        group_by = group_by,
        save_dir = paste0(save_dir, "/", group_by, "/"))
    
    # Run GSEA pipeline
    gsea <- run_gsea_pip(
        res = res,
        org = org,
        order = order,
        save_dir = paste0(save_dir, "/", group_by, "/"))
    
    message("Complete; returning DESeq2 object...")
    return(dds)
}