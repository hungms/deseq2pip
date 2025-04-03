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
#' \dontrun{
#' # Run QC pipeline with default settings
#' dds <- run_qc_pip(dds, experiment = "my_experiment")
#' 
#' # Run QC pipeline with custom settings
#' dds <- run_qc_pip(dds, experiment = "my_experiment",
#'                   org = "mouse", remove_xy = TRUE,
#'                   save_dir_name = "custom_results")
#' }
#' @export
run_qc_pip <- function(
    dds,
    experiment,
    save_dir = getwd(), 
    org = "human",
    remove_xy = FALSE,
    remove_mt = FALSE,
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
    dds <- remove_low_expression(dds, experiment = experiment, quantile = quantile, group_by = group_by, save_dir = save_dir, save_dir_name = save_dir_name)
    check_library(dds, experiment = experiment, save_dir = save_dir, save_dir_name = save_dir_name)

    # Save expression data
    save_expression(dds, experiment = experiment, group_by = group_by, save_dir = save_dir, save_dir_name = save_dir_name)
    
    return(dds)
}

#' Run Sample Distance Pipeline
#'
#' This function runs the distance analysis portion of the RNA-seq analysis pipeline,
#' including generating PCA and distance plots, and saving distance data.
#'
#' @param dds DESeq2 object containing the expression data
#' @param experiment Name of the experiment
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run distance pipeline with default settings
#' dds <- run_dist_pip(dds, experiment = "my_experiment")
#' 
#' # Run distance pipeline with custom settings and batch correction
#' dds <- run_dist_pip(dds, experiment = "my_experiment", 
#'                    batch = "Batch", save_dir_name = "qc_results")
#' }
#' @export
run_dist_pip <- function(
    dds,
    experiment,
    save_dir = getwd(), 
    group_by = "Group1",
    batch = NULL,
    save_dir_name = "qc_results") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(length(experiment) == 1 & is.character(experiment))
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    if(length(batch) > 0){
        stopifnot(batch %in% colnames(colData(dds)))}
    setwd(save_dir)

    # batch correction
    message("Running batch correction...")
    vsd_nobatch <- vst(dds, blind=FALSE)
    colnames(vsd_nobatch) <- colnames(dds)
    
    # Initialize list with non-batch corrected data
    vsd_list <- list(nobatch = vsd_nobatch)
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

        if(method_name != "nobatch"){
            save_dir_name <- "batch_corrected"}

        # Run PCA for additional nfcore groups if present
        nfcore.groups <- colnames(colData(select.vsd))[which(str_detect(colnames(colData(select.vsd)), "^Group"))]
        if (length(nfcore.groups) > 0) {
            for (j in seq_along(nfcore.groups)) {
                run_pca(
                    vsd = select.vsd, 
                    experiment = method_name, 
                    group_by = nfcore.groups[j], 
                    size = 4, 
                    save_dir = save_dir, 
                    save_dir_name = save_dir_name)
            }
        }
            
        # Run PCA and distance analysis
        run_pca(
            vsd = select.vsd, 
            experiment = names(vsd_list)[i], 
            group_by = group_by, 
            size = 4, 
            save_dir = save_dir, 
            save_dir_name = save_dir_name)
        run_distance(
            vsd = select.vsd, 
            experiment = names(vsd_list)[i], 
            save_data = T, 
            save_dir = save_dir, 
            save_dir_name = save_dir_name)
        
        # save vsd object
        save_tsv(assay(select.vsd), experiment = names(vsd_list)[i], tsv_name = "exprs.tsv", save_dir = paste0(save_dir, "/", save_dir_name))
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
    group_by = "Group1",
    save_dir = getwd()) {
    
    # Input validation
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(as.character(design(dds)) != "~ 1")

    # Run differential expression analysis
    res <- run_diffexp_list(dds, org = org, group_by = group_by, save_dir = save_dir)

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
#' @param order Column name to use for ranking genes. Default is "rank"
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
    gsea.df <- read_gsea_tsv_list(merge = TRUE, data_dir = save_dir)
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
#' @param assaytype Type of assay, either "RNA" or "ATAC". Default is "RNA"
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
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
#' \dontrun{
#' # Run complete pipeline with default settings
#' dds <- run_deseq2_pip(dds, experiment = "my_experiment")
#' 
#' # Run pipeline with custom settings
#' dds <- run_deseq2_pip(dds, experiment = "my_experiment",
#'                     org = "mouse", remove_xy = TRUE,
#'                     group_by = "Treatment",
#'                     save_dir_name = "custom_results")
#'
#' # Run pipeline with batch correction
#' dds <- run_deseq2_pip(dds, experiment = "my_experiment",
#'                     batch = "Batch")
#' }
#' @export
run_deseq2_pip <- function(
    dds,
    experiment,
    batch = NULL,
    save_dir = getwd(), 
    org = "human",
    remove_xy = FALSE,
    remove_mt = FALSE,
    group_by = "Group1",
    quantile = 0.05,
    order = "rank",
    save_dir_name = "qc_results",
    assaytype = "RNA") {
    
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
    stopifnot(assaytype %in% c("RNA", "ATAC"))
    if(length(batch) > 0){
        stopifnot(batch %in% colnames(colData(dds)))}
    if(assaytype == "ATAC"){
        stopifnot("TSS" %in% colnames(rowData(dds)))}


    # Run quality control pipeline
    message("Running quality control pipeline...")
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

    # Run distance pipeline
    message("Running distance analysis pipeline...")
    dds <- run_dist_pip(
        dds = dds, 
        experiment = experiment, 
        batch = batch,
        save_dir = save_dir, 
        group_by = group_by, 
        save_dir_name = save_dir_name)
    
    # Run differential expression pipeline
    message("Running differential expression analysis...")
    res <- run_diffexp_pip(
        dds = dds,
        org = org,
        group_by = group_by,
        save_dir = paste0(save_dir, "/", group_by, "/"))
    
    # Run GSEA pipeline for RNA-seq
    if(assaytype == "RNA"){
        message("Running GSEA pipeline...")
        gsea <- run_gsea_pip(
            res = res,
            org = org,
            order = order,
            save_dir = paste0(save_dir, "/", group_by, "/"))

        message("Complete; returning DESeq2 object...")}

    # Run DESeq2 and GSEA pipeline for ATAC-seq TSS peaks
    if(assaytype == "ATAC"){
        message("Repeating analysis for TSS peaks...")

        # Subset dds to TSS peaks only
        dds.tss <- getTSS(dds)
        dds.tss[[paste0(group_by, "_TSS")]] <- factor(dds.tss[[group_by]], levels = levels(dds.tss[[group_by]]))
        design(dds.tss) <- as.formula(paste0("~", group_by, "_TSS"))

        # Run quality control pipeline
        dds.tss <- run_dist_pip(
            dds = dds.tss, 
            experiment = experiment, 
            batch = batch,
            save_dir = save_dir, 
            group_by = paste0(group_by, "_TSS"), 
            save_dir_name = paste0(save_dir_name, "_TSS"))
        
        # Run differential expression pipeline
        res.tss <- run_diffexp_pip(
            dds = dds.tss,
            org = org,
            group_by = paste0(group_by, "_TSS"),
            save_dir = paste0(save_dir, "/", group_by, "_TSS/"))
        
        # Run GSEA pipeline
        gsea.tss <- run_gsea_pip(
            res = res.tss,
            org = org,
            order = order,
            save_dir = paste0(save_dir, "/", group_by, "_TSS/"))

        # Repeat for PCA and distance analysis for TSS peaks only
        message("TSS peak analysis complete; returning DESeq2 object...")
    }
    
    return(dds)
}
