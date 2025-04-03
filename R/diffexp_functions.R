#' Run Differential Expression Analysis
#'
#' This function performs differential expression analysis using DESeq2 for a comparison
#' between two groups. It calculates fold changes, p-values, and adjusted p-values for
#' each gene, and ranks genes based on their statistical and biological significance.
#'
#' @param dds DESeq2 object containing the expression data
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param save_data Logical. If TRUE, saves results to TSV file. Default is TRUE
#' @param save_dir Directory to save the results. Default is current working directory
#' @param ... Additional arguments passed to DESeq2::DESeq()
#' @return A data frame containing differential expression results with columns:
#'         - gene: gene identifier
#'         - comparison: name of the comparison
#'         - log2FoldChange: log2 fold change between groups
#'         - padj: adjusted p-value
#'         - rank: combined score (-log10(padj) * log2FoldChange)
#' @examples
#' \dontrun{
#' # Basic differential expression analysis
#' res <- run_diffexp(dds)
#' 
#' # Differential expression analysis with custom parameters
#' res <- run_diffexp(dds, group_by = "Treatment", save_data = TRUE)
#' }
#' @export
run_diffexp <- function(dds, org = "human", group_by = "Group1", save_data = T, save_dir = getwd(), ...){
    # make sure there are 2 levels
    meta <- colData(dds) %>% as.data.frame(.)
    row.meta <- rowData(dds) %>% as.data.frame(.)
    stopifnot(is.factor(meta[[group_by]]))
    stopifnot(length(levels(meta[[group_by]])) == 2)

    # set comparisons
    group.lv <- levels(meta[[group_by]])
    comparison <- paste0(group.lv[2], "_vs_", group.lv[1])

    # run deseq2
    dds <- DESeq(dds, ...)
    res <- results(dds)

    # shrink lfc by ashr
    res <- lfcShrink(dds, res = res, type = "ashr") %>% 
        as.data.frame(.)
    
    # add row metadata
    row.meta[which(colnames(row.meta) %in% colnames(res))] <- NULL

    if(ncol(row.meta) > 0){
        message("Combining row metadata to results...")
        if(all(rownames(res) == rownames(assay(dds)))){
            res <- res[rownames(assay(dds)),]
            res <- cbind(res, row.meta)}
    }

    # modify columns
    res <- res %>%
        mutate(
            comparison = comparison,
            rank = -log10(padj)*log2FoldChange,
            rank = ifelse(is.infinite(rank), 0, rank)
            ) %>%
        arrange(desc(rank))
    
    if(!"gene" %in% colnames(res)){
        res <- res %>%
            mutate(gene = rownames(.))}
    
    res <- run_annotation(res, org = org, gene_column = "gene")

    # save data
    if(save_data){
        save_tsv(res, tsv_name = "diffexp_deseq2_wald.tsv", save_dir = paste0(save_dir, comparison, sep = "/"))}

    return(res)
}


#' Run Differential Expression Analysis for All Comparisons
#'
#' This function performs differential expression analysis for all possible pairwise
#' comparisons between groups in the experiment.
#'
#' @param dds DESeq2 object containing the expression data
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param ... Additional arguments passed to run_diffexp()
#' @return A merged data frame containing differential expression results for all comparisons
#' @examples
#' \dontrun{
#' # Run differential expression for all comparisons
#' res <- run_diffexp_list(dds)
#' 
#' # Run with custom parameters
#' res <- run_diffexp_list(dds, group_by = "Treatment")
#' }
#' @export
run_diffexp_list <- function(dds,  group_by = "Group1", ...){
    message("Running DESeq2 differential expression...")
    stopifnot(is.factor(colData(dds)[[group_by]]))
    meta <- colData(dds) %>% as.data.frame(.)
    group.lv <- levels(meta[[group_by]])

    res.list <- list()
    n_comparisons <- choose(length(group.lv), 2)
    pb <- txtProgressBar(min = 0, max = n_comparisons, style = 3)
    counter <- 0

    for(g1 in 1:length(group.lv)){
        for(g2 in 1:length(group.lv)){
            if(g1 >= g2){next}
            
            samples <- meta %>%
                filter(str_detect(!!sym(group_by),  paste0("^", group.lv[g1], "$|^", group.lv[g2], "$"))) %>%
                rownames(.)

            dds.selected <- dds[,samples]
            colData(dds.selected)[[group_by]] <- factor(colData(dds.selected)[[group_by]], levels = c(group.lv[g1], group.lv[g2]))
            comparison <- paste0(group.lv[g2], "_vs_", group.lv[g1])

            res.list[[length(res.list) + 1]] <- run_diffexp(dds.selected, group_by = group_by, ...)
            counter <- counter + 1
            setTxtProgressBar(pb, counter)
        }
    }
    close(pb)
    res.list <- bind_rows(res.list)
    return(res.list)
}

#' Read Differential Expression Results
#'
#' This function reads differential expression results from multiple comparison directories
#' and optionally merges them into a single data frame.
#'
#' @param data_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param merge Logical. If TRUE, returns a single merged data frame. If FALSE, returns a list of data frames
#' @return Either a merged data frame or a list of data frames containing differential expression results
#' @examples
#' \dontrun{
#' # Read and merge all results
#' merged_results <- read_diffexp_list(data_dir = "results", merge = TRUE)
#' 
#' # Read results as a list
#' results_list <- read_diffexp_list(data_dir = "results", merge = FALSE)
#' }
#' @export
read_diffexp_list <- function(data_dir = getwd(), merge = T){

    comparison <- list.files(data_dir)
    files <- c()
    for(i in seq_along(comparison)){
        files <- c(files, paste0(comparison[i], "/", list.files(paste0(data_dir, "/", comparison[i]))))}
    files <- files[which(str_detect(files, "/diffexp_deseq2_wald.tsv"))]
    res.list <- list()

    for(i in seq_along(files)){
        res.list[[i]] <- read.table(paste0(data_dir, "/", files[i]), sep = "\t", header = T)

        if(!"comparison" %in% colnames(res.list[[i]])){
            res.list[[i]]$comparison <- files[i]
            names(res.list)[i] <- files[i]}
        else{
            names(res.list)[i] <- unique(res.list[[i]]$comparison)}}

    if(merge){
        res.list <- bind_rows(res.list)}

    return(res.list)
}


#' Generate MA Plots for All Comparisons
#'
#' This function creates MA plots for all comparisons in the differential
#' expression results.
#'
#' @param res Differential expression result data frame from run_diffexp()
plot_ma_list <- function(res, ...){
    message("Generating MA plots...")
    stopifnot(is.data.frame(res))
    stopifnot(is.data.frame(res) & all(c("baseMean", "log2FoldChange", "padj", "gene", "comparison") %in% colnames(res)))
    res.list <- split(res, res$comparison)

    ma.plot.list <- list()
    pb <- txtProgressBar(min = 0, max = length(res.list), style = 3)
    
    for(i in seq_along(res.list)){
        ma.plot <- plot_ma(res.list[[i]], ...)
        ma.plot.list <- c(ma.plot.list, ma.plot)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    return(ma.plot.list)
}


#' Generate Volcano Plots for All Comparisons
#'
#' This function creates volcano plots for all comparisons in the differential
#' expression results.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param ... Additional arguments passed to plot_volcano()
#' @return A list of volcano plots for each comparison
#' @examples
#' \dontrun{
#' # Generate volcano plots for all comparisons
#' plots <- plot_volcano_list(res)
#' 
#' # Generate with custom parameters
#' plots <- plot_volcano_list(res, n = 30, fc.thresh = 2)
#' }
#' @export
plot_volcano_list <- function(res, ...){
    message("Generating volcano plots...")
    stopifnot(is.data.frame(res))
    stopifnot(c("padj", "gene", "log2FoldChange", "comparison") %in% colnames(res))
    res.list <- split(res, res$comparison)

    volcano.plot.list <- list()
    pb <- txtProgressBar(min = 0, max = length(res.list), style = 3)
    
    for(i in seq_along(res.list)){
        volcano.plot <- plot_volcano(res.list[[i]], ...)
        volcano.plot.list <- c(volcano.plot.list, volcano.plot)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    return(volcano.plot.list)
}