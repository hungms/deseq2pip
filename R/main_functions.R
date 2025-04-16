#' Run Complete RNA-seq Pipeline
#'
#' This function runs a complete RNA-seq analysis pipeline including quality control,
#' differential expression analysis, and gene set enrichment analysis for all possible
#' comparisons in the experiment. It generates various plots and saves results in
#' organized directories.
#'
#' @param dds DESeq2 object containing the expression data
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param quantile Quantile threshold for filtering lowly expressed genes. Default is 0.05
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run complete pipeline with default settings
#' dds <- run_rna_pip(dds)
#' 
#' # Run pipeline with custom settings
#' dds <- run_rna_pip(dds, org = "mouse", remove_xy = TRUE,
#'                     group_by = "Treatment",
#'                     save_dir_name = "custom_results")
#'
#' # Run pipeline with batch correction
#' dds <- run_rna_pip(dds, batch = "Batch")
#' }
#' @export
run_rna_pip <- function(
    dds,
    group_by = "Group1",
    pals = NULL,
    batch = NULL,
    save_dir = getwd(), 
    org = "human",
    remove_xy = FALSE,
    remove_mt = FALSE,
    quantile = 0.05,
    order = "pxfc",
    save_dir_name = "qc_results") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(org %in% c("mouse", "human"))
    stopifnot(is.logical(remove_xy))
    stopifnot(is.logical(remove_mt))
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(is.numeric(quantile) & quantile < 0.2)
    stopifnot(order %in% c("log2FoldChange", "padj", "pxfc"))
    stopifnot(as.character(design(dds)) != "~ 1")
    if(length(batch) > 0){
        stopifnot(batch %in% colnames(colData(dds)))}
    if(length(pals) > 0){
        stopifnot(all(colData(dds)[[group_by]] %in% names(pals)))}

    # Display package version
    message("Running RNA-seq pipeline with DESeq2pip v", packageVersion("deseq2pip"))

    # Run quality control pipeline
    message("Running quality control pipeline...")
    dds <- run_qc_pip(
        dds = dds,
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
        pals = pals,
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
        order = order,
        save_dir = paste0(save_dir, "/", group_by, "/"))
    
    # Run GSEA pipeline for RNA-seq
    message("Running GSEA pipeline...")
        gsea <- run_gsea_pip(
        res = res,
        org = org,
        save_dir = paste0(save_dir, "/", group_by, "/"))

    message("RNAseq pipeline complete; returning DESeq2 object...")
    return(dds)
}



#' Run Complete ATAC-seq Pipeline
#'
#' This function runs a complete ATAC-seq analysis pipeline including quality control,
#' differential expression analysis, and gene set enrichment analysis for all possible
#' comparisons in the experiment. It generates various plots and saves results in
#' organized directories.
#'
#' @param dds DESeq2 object containing the expression data
#' @param assaytype Type of assay, either "RNA" or "ATAC". Default is "RNA"
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
#' @param repeat_for_TSS Logical. If TRUE, repeats the analysis for TSS peaks. Default is TRUE
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param quantile Quantile threshold for filtering lowly expressed genes. Default is 0.05
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param save_dir_name Name of the subdirectory to save files in. Default is "qc_results"
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run complete pipeline with default settings
#' dds <- run_atac_pip(dds)
#' 
#' # Run pipeline with custom settings
#' dds <- run_atac_pip(dds, org = "mouse", remove_xy = TRUE,
#'                     group_by = "Treatment",
#'                     save_dir_name = "custom_results")
#'
#' # Run pipeline with batch correction
#' dds <- run_atac_pip(dds, batch = "Batch")
#' }
#' @export
run_atac_pip <- function(
    dds,
    group_by = "Group1",
    repeat_for_TSS = TRUE,
    pals = NULL,
    batch = NULL,
    save_dir = getwd(), 
    org = "human",
    remove_xy = FALSE,
    remove_mt = FALSE,
    quantile = 0.05,
    order = "pxfc",
    save_dir_name = "qc_results",
    assaytype = "RNA") {
    
    # Input validation
    stopifnot(dir.exists(save_dir))
    stopifnot(all(c("peaks", "gene", "Annotation", "TSS") %in% colnames(rowData(dds))))
    stopifnot(org %in% c("mouse", "human"))
    stopifnot(is.logical(remove_xy))
    stopifnot(is.logical(remove_mt))
    stopifnot(group_by %in% colnames(colData(dds)))
    stopifnot(is.factor(colData(dds)[[group_by]]))
    stopifnot(is.numeric(quantile) & quantile < 0.2)
    stopifnot(order %in% c("log2FoldChange", "padj", "pxfc"))
    stopifnot(as.character(design(dds)) != "~ 1")
    if(length(batch) > 0){
        stopifnot(batch %in% colnames(colData(dds)))}
    if(length(pals) > 0){
        stopifnot(all(colData(dds)[[group_by]] %in% names(pals)))}

    # Display package version
    message("Running ATAC-seq pipeline version with DESeq2pip v", packageVersion("deseq2pip"))

    # Run quality control pipeline
    message("Running quality control pipeline...")
    dds <- run_qc_pip(
        dds = dds,
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
        pals = pals,
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
        order = order,
        save_dir = paste0(save_dir, "/", group_by, "/"))

    # Generate ATACseq annotation plots for all comparisons
    message("Generating ATACseq annotation plots...")
    plist <- plot_atac_annot_list(dds, res, save = TRUE, save_dir = paste0(save_dir, "/", group_by))

    if(repeat_for_TSS){
        # Run DESeq2 and GSEA pipeline for ATAC-seq TSS peaks
        message("Repeating analysis for TSS peaks...")
        dds.tss <- getTSS(dds)
        dds.tss[[paste0(group_by, "_TSS")]] <- factor(dds.tss[[group_by]], levels = levels(dds.tss[[group_by]]))
        design(dds.tss) <- as.formula(paste0("~", group_by, "_TSS"))

        # Run quality control pipeline
        message("Running distance analysis pipeline for TSS peaks...")
        dds.tss <- run_dist_pip(
            dds = dds.tss, 
            batch = batch,
            save_dir = save_dir, 
            group_by = paste0(group_by, "_TSS"), 
            save_dir_name = paste0(save_dir_name, "_TSS"))
            
        # Run differential expression pipeline
        message("Running differential expression analysis for TSS peaks...")
        res.tss <- run_diffexp_pip(
            dds = dds.tss,
            org = org,
            order = order,
            group_by = paste0(group_by, "_TSS"),
            save_dir = paste0(save_dir, "/", group_by, "_TSS/"))
            
        # Run GSEA pipeline
        message("Running GSEA pipeline for TSS peaks...")
        gsea.tss <- run_gsea_pip(
            res = res.tss,
            org = org,
            save_dir = paste0(save_dir, "/", group_by, "_TSS/"))}

    # Repeat for PCA and distance analysis for TSS peaks only
    message("ATACseq analysis complete; returning DESeq2 object...")
    
    return(dds)
}
