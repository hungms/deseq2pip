#' Run Complete RNA-seq Pipeline
#'
#' This function runs a complete RNA-seq analysis pipeline including quality control,
#' differential expression analysis, and gene set enrichment analysis for all possible
#' comparisons in the experiment. It generates various plots and saves results in
#' organized directories.
#'
#' @param dds DESeq2 object
#' @param org The organism to use, either "human" or "mouse".
#' @param var Column name in colData(dds) to use for grouping.
#' @param design Design formula for the DESeq2 object.
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL
#' @param comparisons Vector of comparisons to run. Default is NULL
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param min_count Minimum count threshold for filtering lowly expressed genes. Default is 10
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run complete pipeline with default settings
#' dds <- run_rna_pip(dds)
#' 
#' # Run pipeline with custom settings
#' dds <- run_rna_pip(dds, 
#'                    org = "mouse", 
#'                    remove_xy = TRUE,
#'                    var = "Treatment",
#'                    save_dir = "output")
#'
#' }
#' @export
#' 
run_rna_pip <- function(
    dds,
    org,
    var,
    design,
    batch = NULL,
    comparisons = NULL,
    pals = NULL,
    remove_xy = FALSE,
    remove_mt = FALSE,
    min_count = 10,
    order = "pxfc",
    save_dir = getwd()) {

    log_output(match.call(), save_dir, expr = quote({
            
        # validations
        message("Running RNA-seq pipeline with DESeq2pip v", packageVersion("deseq2pip"))
        dds <- validate_var(dds, var)
        dds <- validate_design(dds, design)
        org <- validate_org(org)
        min_count <- validate_min_count(min_count)
        order <- validate_order(order)
        validate_logical(c(remove_xy, remove_mt))
        if(!is.null(batch)) dds <- validate_batch(dds, batch)
        if(!is.null(pals)) pals <- validate_pals(dds, var, pals)
        comparisons <- validate_comparisons(dds, var, comparisons)

        # Run quality control pipeline
        dds <- run_qc_pip(
            dds = dds,
            org = org,
            var = var,
            remove_xy = remove_xy,
            remove_mt = remove_mt,
            min_count = min_count,
            save_dir = save_dir)

        # Run distance pipeline
        dds <- run_dist_pip(
            dds = dds, 
            var = var, 
            batch = batch,
            pals = pals,
            save_dir = save_dir)

        # Run differential expression analysis for all comparisons
        res <- run_diffexp_pip(
            dds = dds,
            org = org,
            var = var,
            design = design,
            comparisons = comparisons,
            order = order,
            save_dir = paste0(save_dir, "/", var))
        
        # Run gene set enrichment analysis for all comparisons
        gsea <- run_gsea_pip(
            res = res,
            org = org,
            save_dir = paste0(save_dir, "/", var))

        # Run enrichment map analysis for all comparisons
        enrichmentmap_pip(dds, org, var, group_dir = paste0(save_dir, "/", var))

        message("RNA-seq pipeline complete; returning DESeq2 object...")
        return(dds)
    }))
}



#' Run Complete ATAC-seq Pipeline
#'
#' This function runs a complete ATAC-seq analysis pipeline including quality control,
#' differential expression analysis, and gene set enrichment analysis for all possible
#' comparisons in the experiment. It generates various plots and saves results in
#' organized directories.
#'
#' @param dds DESeq2 object
#' @param org The organism to use, either "human" or "mouse".
#' @param var Column name in colData(dds) to use for grouping.
#' @param design Design formula for the DESeq2 object.
#' @param comparisons Vector of comparisons to run. Default is NULL
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param min_count Minimum count threshold for filtering lowly expressed genes. Default is 10
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param TSS Logical. If TRUE, repeats the analysis for TSS peaks. Default is TRUE  
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return The processed DESeq2 object
#' @examples
#' \dontrun{
#' # Run complete pipeline with default settings
#' dds <- run_atac_pip(dds)
#' 
#' # Run pipeline with custom settings
#' dds <- run_atac_pip(dds, 
#'                     org = "mouse", 
#'                     remove_xy = TRUE,
#'                     var = "Treatment",
#'                     save_dir = "output")
#'
#' }
#' @export
run_atac_pip <- function(
    dds,
    org,
    var,
    design,
    batch = NULL,
    comparisons = NULL,
    pals = NULL,
    remove_xy = FALSE,
    remove_mt = FALSE,
    min_count = 10,
    order = "pxfc",
    TSS = TRUE,
    save_dir = getwd()) {

    log_output(match.call(), save_dir, expr = quote({
        
        # validations
        message("Running ATAC-seq pipeline with DESeq2pip v", packageVersion("deseq2pip"))
        dds <- validate_dds_atac(dds)
        dds <- validate_var(dds, var)
        dds <- validate_design(dds, design)
        org <- validate_org(org)
        min_count <- validate_min_count(min_count)
        order <- validate_order(order)
        validate_logical(c(remove_xy, remove_mt, TSS))
        if(!is.null(batch)) dds <- validate_batch(dds, batch)
        if(!is.null(pals)) pals <- validate_pals(dds, var, pals)
        comparisons <- validate_comparisons(dds, var, comparisons)

        # Run quality control pipeline
        dds <- run_qc_pip(
            dds = dds,
            org = org,
            var = var,
            remove_xy = remove_xy,
            remove_mt = remove_mt,
            min_count = min_count,
            save_dir = save_dir)

        # Run distance pipeline
        dds <- run_dist_pip(
            dds = dds, 
            var = var, 
            batch = batch,
            pals = pals,
            save_dir = save_dir)

        # Run differential expression analysis for all comparisons  
        res <- run_diffexp_pip(
            dds = dds,
            org = org,
            var = var,
            design = design,
            comparisons = comparisons,
            order = order,
            save_dir = paste0(save_dir, "/", var))

        # Generate peak annotation plots for all comparisons
        plot_peak_annot_pip(dds, res, save_dir = paste0(save_dir, "/", var))

        if(TSS){
            # Get TSS peaks
            tss_dds <- getTSS(dds)
            tss_var <- paste0(var, "_TSS")
            tss_dds[[tss_var]] <- factor(tss_dds[[var]], levels = levels(tss_dds[[var]]))
            tss_design <- gsub(var, tss_var, design)
            tss_dds <- validate_design(tss_dds, tss_design)
            design(tss_dds) <- as.formula(tss_design)


            # Run differential expression pipeline for TSS peaks
            tss_res <- run_diffexp_pip(
                dds = tss_dds,
                org = org,
                var = tss_var,
                design = tss_design,
                comparisons = comparisons,
                order = order,
                save_dir = paste0(save_dir, "/", var, "_TSS"))
            
            # Run gene set enrichment pipeline for TSS peaks
            tss_gsea <- run_gsea_pip(
                res = tss_res,
                org = org,
                save_dir = paste0(save_dir, "/", var, "_TSS"))

            # Run enrichment map pipeline for TSS peaks
            enrichmentmap_pip(tss_dds, org, tss_var, group_dir = paste0(save_dir, "/", var, "_TSS"))
        }

        message("ATAC-seq pipeline complete; returning DESeq2 object...")
        return(dds)
        }))
}