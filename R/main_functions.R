#' Run Complete RNA-seq Pipeline
#'
#' This function runs a complete RNA-seq analysis pipeline including quality control,
#' differential expression analysis, and gene set enrichment analysis for all possible
#' comparisons in the experiment. It generates various plots and saves results in
#' organized directories.
#'
#' @param dds DESeq2 object
#' @param org The organism to use, either "human" or "mouse".
#' @param group_by Column name in colData(dds) to use for grouping.
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param min_count Minimum count threshold for filtering lowly expressed genes. Default is 10
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL 
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
#'                    group_by = "Treatment",
#'                    save_dir = "output")
#'
#' # Run pipeline with batch correction
#' dds <- run_rna_pip(dds, batch = "Batch")
#' }
#' @export
run_rna_pip <- function(
    dds,
    org,
    group_by,
    remove_xy = FALSE,
    remove_mt = FALSE,
    min_count = 10,
    pals = NULL,
    batch = NULL,
    order = "pxfc",
    save_dir = getwd()) {

    log_output(match.call(), save_dir, expr = quote({
            
        # validations
        validate_logical(c(remove_xy, remove_mt))
        dir.create(save_dir, recursive = TRUE)
        validate_paths(save_dir)
        dds <- validate_dds_group_by(dds, group_by)
        org <- validate_org(org)
        min_count <- validate_min_count(min_count)
        pals <- validate_pals(dds, group_by, pals)
        dds <- validate_batch(dds, batch)
        order <- validate_order(order)

        message("Running RNA-seq pipeline with DESeq2pip v", packageVersion("deseq2pip"))

        # Run quality control pipeline
        message("Running quality control pipeline...")
        dds <- run_qc_pip(
            dds = dds,
            org = org,
            remove_xy = remove_xy,
            remove_mt = remove_mt,
            group_by = group_by,
            min_count = min_count,
            save_dir = save_dir)

        # Run distance pipeline
        message("Running distance analysis pipeline...")
        dds <- run_dist_pip(
            dds = dds, 
            group_by = group_by, 
            pals = pals,
            batch = batch,
            save_dir = save_dir)

        # Run pairwise differential expression pipeline
        message("Running PAIRWISE differential expression analysis...")
        res <- run_diffexp_pip(
            dds = dds,
            org = org,
            group_by = group_by,
            order = order,
            one_to_all = FALSE,
            save_dir = save_dir)
        
        # Run pairwise GSEA pipeline
        message("Running PAIRWISE gene set enrichment pipeline...")
        group_save_dir <- paste0(save_dir, "/pairwise_", group_by, "/")
        gsea <- run_gsea_pip(
            res = res,
            org = org,
            group_save_dir = group_save_dir)

        # Run PAIRWISE enrichment map pipeline
        message("Running PAIRWISE enrichment map pipeline...")
        enrichmentmap_pip(dds, org, group_by, group_dir = group_save_dir)

        # Whether to run ONE-TO-ALL comparisons
        if(length(unique(colData(dds)[[group_by]])) > 2){

            # Run ONE-TO-ALL differential expression pipeline
            message("Running ONE-TO-ALL differential expression analysis...")
            res <- run_diffexp_pip(
                dds = dds,
                org = org,
                group_by = group_by,
                order = order,
                one_to_all = TRUE,
                save_dir = save_dir)
            
            # Run ONE-TO-ALL GSEA pipeline
            message("Running ONE-TO-ALL gene set enrichment pipeline...")
            group_save_dir <- paste0(save_dir, "/one-to-all_", group_by, "/")
            gsea <- run_gsea_pip(
                res = res,
                org = org,
                group_save_dir = group_save_dir)

            # Run ONE-TO-ALL enrichment map pipeline
            message("Running ONE-TO-ALL enrichment map pipeline...")
            enrichmentmap_pip(dds, org, group_by, group_dir = group_save_dir)}

        message("RNAseq pipeline complete; returning DESeq2 object...")
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
#' @param group_by Column name in colData(dds) to use for grouping.
#' @param remove_xy Logical. If TRUE, removes genes on X and Y chromosomes. Default is FALSE
#' @param remove_mt Logical. If TRUE, removes mitochondrial genes. Default is FALSE
#' @param min_count Minimum count threshold for filtering lowly expressed genes. Default is 10
#' @param pals Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param batch Column name in colData(dds) to use for batch correction. Default is NULL 
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
#'                     group_by = "Treatment",
#'                     save_dir = "output")
#'
#' # Run pipeline with batch correction
#' dds <- run_atac_pip(dds, batch = "Batch")
#' }
#' @export
run_atac_pip <- function(
    dds,
    org,
    group_by,
    remove_xy = FALSE,
    remove_mt = FALSE,
    min_count = 10,
    pals = NULL,
    batch = NULL,
    order = "pxfc",
    TSS = TRUE,
    save_dir = getwd()) {

    log_output(match.call(), save_dir, expr = quote({
        
        # validations
        validate_logical(c(remove_xy, remove_mt))
        dir.create(save_dir, recursive = TRUE)
        validate_paths(save_dir)
        dds <- validate_dds_atac(dds)
        dds <- validate_dds_group_by(dds, group_by)
        org <- validate_org(org)
        min_count <- validate_min_count(min_count)
        pals <- validate_pals(dds, group_by, pals)
        dds <- validate_batch(dds, batch)
        order <- validate_order(order)

        # Display package version
        message("Running ATAC-seq pipeline with DESeq2pip v", packageVersion("deseq2pip"))

        # Run quality control pipeline
        message("Running quality control pipeline...")
        dds <- run_qc_pip(
            dds = dds,
            org = org,
            remove_xy = remove_xy,
            remove_mt = remove_mt,
            group_by = group_by,
            min_count = min_count,
            save_dir = save_dir)

        # Run distance pipeline
        message("Running distance analysis pipeline...")
        dds <- run_dist_pip(
            dds = dds, 
            group_by = group_by, 
            pals = pals,
            batch = batch,
            save_dir = save_dir)

        # Run pairwise differential expression pipeline
        message("Running PAIRWISE differential expression analysis for ALL PEAKS...")
        res <- run_diffexp_pip(
            dds = dds,
            org = org,
            group_by = group_by,
            order = order,
            one_to_all = FALSE,
            save_dir = save_dir)

        # Generate one-to-all ATACseq annotation plots for all comparisons for all peaks
        message("Generating PAIRWISE peak annotation plots...")
        group_save_dir <- paste0(save_dir, "/pairwise_", group_by, "/")
        plot_peak_annot_pip(dds, res, group_save_dir = group_save_dir)

        # Whether to run ONE-TO-ALL comparisons
        if(length(unique(colData(dds)[[group_by]])) > 2){

            # Run ONE-TO-ALL differential expression pipeline for ALL PEAKS
            message("Running ONE-TO-ALL differential expression analysis for ALL PEAKS...")
            res <- run_diffexp_pip(
                dds = dds,
                org = org,
                group_by = group_by,
                order = order,
                one_to_all = TRUE,
                save_dir = save_dir)

            # Generate one-to-all ATACseq annotation plots for all comparisons for ALL PEAKS
            message("Generating ONE-TO-ALL peak annotation plots for ALL PEAKS...")
            group_save_dir <- paste0(save_dir, "/one-to-all_", group_by, "/")
            plot_peak_annot_pip(dds, res, group_save_dir = group_save_dir)}

        if(TSS){
            tss.dds <- getTSS(dds)
            tss_group_by <- paste0("TSS_", group_by)
            tss.dds[[tss_group_by]] <- factor(tss.dds[[group_by]], levels = levels(tss.dds[[group_by]]))
            des <- paste0(as.character(design(dds)), collapse = " ")
            des <- gsub(group_by, tss_group_by, des)
            design(tss.dds) <- as.formula(des)

            # Run PAIRWISE differential expression pipeline for TSS PEAKS
            message("Running PAIRWISE differential expression analysis for TSS PEAKS...")
            res <- run_diffexp_pip(
                dds = tss.dds,
                org = org,
                group_by = tss_group_by,
                order = order,
                one_to_all = FALSE,
                save_dir = save_dir)
            
            # Run PAIRWISE GSEA pipeline for TSS PEAKS
            message("Running PAIRWISE gene set enrichment pipeline for TSS PEAKS...")
            group_save_dir <- paste0(save_dir, "/pairwise_", tss_group_by, "/")
                gsea <- run_gsea_pip(
                    res = res,
                    org = org,
                    group_save_dir = group_save_dir)

            # Run PAIRWISE enrichment map pipeline for TSS PEAKS
            message("Running PAIRWISE enrichment map pipeline for TSS PEAKS...")
            enrichmentmap_pip(tss.dds, org, tss_group_by, group_dir = group_save_dir)


            # Whether to run ONE-TO-ALL comparisons for TSS PEAKS
            if(length(unique(colData(dds)[[group_by]])) > 2){

                # Run ONE-TO-ALL differential expression pipeline for TSS PEAKS
                message("Running ONE-TO-ALL differential expression analysis for TSS PEAKS...")
                res <- run_diffexp_pip(
                    dds = tss.dds,
                    org = org,
                    group_by = tss_group_by,
                    order = order,
                    one_to_all = TRUE,
                    save_dir = save_dir)
                
                # Run ONE-TO-ALL GSEA pipeline for TSS PEAKS
                message("Running ONE-TO-ALL gene set enrichment pipeline for TSS PEAKS...")
                group_save_dir <- paste0(save_dir, "/one-to-all_", tss_group_by, "/")
                gsea <- run_gsea_pip(
                    res = res,
                    org = org,
                    group_save_dir = group_save_dir)

                # Run ONE-TO-ALL enrichment map pipeline for TSS PEAKS
                message("Running ONE-TO-ALL enrichment map pipeline for TSS PEAKS...")
                enrichmentmap_pip(tss.dds, org, tss_group_by, group_dir = group_save_dir)}
        }

        message("ATACseq pipeline complete; returning DESeq2 object...")
        return(dds)
        }))
}