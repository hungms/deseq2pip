#' Validate DESeq2 object
#'
#' This function validates a DESeq2 object.
#'
#' @param dds A DESeq2 object.
#' @return A DESeq2 object.
#' @export
validate_dds <- function(dds){
    if(!"gene" %in% colnames(rowData(dds))){
        message("No gene column detected in rowData, adding rownames as gene column")
        rowData(dds)$gene <- rownames(dds)}
    if(class(dds) != "DESeqTransform"){
        stopifnot(as.character(design(dds)) != "~ 1")}
    return(dds)}

#' Validate DESeq2 object for ATAC-seq
#'
#' This function validates a DESeq2 object for ATAC-seq.
#'
#' @param dds A DESeq2 object.
#' @return A DESeq2 object.
#' @export
validate_dds_atac <- function(dds){
    if(!"peaks" %in% colnames(rowData(dds))){
        if(all(str_detect(rownames(dds), "Intervals"))){
            message("No peaks column detected in rowData, adding rownames as peaks column")
            rowData(dds)$peaks <- rownames(dds)}
        else{
            stop("Please provide a peaks column in rowData manually")}}
    else{
        stopifnot("Please make sure that the rowData contains peaks, TSS, and Annotation columns" = all(c("peaks", "TSS", "Annotation") %in% colnames(rowData(dds))))}
    validate_logical(rowData(dds)$TSS)
    return(dds)}

#' Validate DESeq2 object for RNA-seq
#'
#' This function validates a DESeq2 object for RNA-seq.
#'
#' @param dds A DESeq2 object.
#' @return A DESeq2 object.
#' @export          
validate_dds_group_by <- function(dds, group_by){

    # check if group_by is a factor
    if(!is.factor(colData(dds)[[group_by]])){
        message("group_by is not a factor, converting to factor")
        colData(dds)[[group_by]] <- factor(colData(dds)[[group_by]], levels = unique(colData(dds)[[group_by]]))}

    if(class(dds) != "DESeqTransform"){
    if(!str_detect(paste0(as.character(design(dds)), collapse = " "), group_by)){
        message("group_by is not in the design, adding to design")
        if(paste0(as.character(design(dds)), collapse = " ") == "~ 1"){
            design(dds) <- as.formula(paste0("~ ", group_by))
        } else {
            design(dds) <- as.formula(paste0(as.character(design(dds)), " + ", group_by))}}
    }

    dds <- validate_dds(dds)
    return(dds)}

#' Validate DESeq2 object for comparison
#'
#' This function validates a DESeq2 object for comparison.
#'
#' @param dds A DESeq2 object.
#' @return A DESeq2 object.
#' @export
validate_dds_comparison <- function(dds, group_by, comparison){
    stopifnot("Please make sure that each group variable is present in group_by" = all(str_split(comparison, "_vs_") %>% unlist(.) %in% colData(dds)[[group_by]]))
    dds <- validate_dds_group_by(dds, group_by)
    return(dds)}

#' Validate methods
#' 
#' This function validates methods.
#' 
#' @param method A character string of the method.  
#' @return A character string of the method.
#' @export
validate_method <- function(method){
    if(!is.null(method)){
        stopifnot("Please make sure that method is either limma or combat" = method %in% c("limma", "combat"))}
    return(method)}

#' Validate differential expression results
#' 
#' This function validates differential expression results.
#'
#' @param res A data frame of differential expression results.
#' @return A data frame of differential expression results.
#' @export
validate_res <- function(res){  
    stopifnot("Please make sure that the differential expression results is a data frame and contains baseMean, log2FoldChange, padj, gene, rank, and comparison columns" = is.data.frame(res) & all(c("baseMean", "log2FoldChange", "padj", "gene", "rank", "comparison") %in% colnames(res)))
    for(c in c("baseMean", "log2FoldChange", "padj", "rank")) {
        message(paste0("Converting ", c, " to numeric"))
        res[[c]] <- as.numeric(res[[c]])}
    return(res)}

#' Validate differential expression results for comparison
#' 
#' This function validates differential expression results for comparison.
#'
#' @param res A data frame of differential expression results.
#' @return A data frame of differential expression results.
#' @export
validate_res_comparison <- function(res){
    res <- validate_res(res)
    stopifnot("Please make sure that the differential expression results contains only one comparison pair" = length(unique(res$comparison)) == 1)
    return(res)}

#' Validate organism
#' 
#' This function validates an organism.
#'
#' @param org A character string of the organism.
#' @return A character string of the organism.
#' @export
validate_org <- function(org){
    stopifnot("Please make sure that the organism is either human or mouse" = org %in% c("human", "mouse"))
    return(org)}

#' Validate logical
#' 
#' This function validates a logical.
#'
#' @param logical A logical.
#' @return A logical.
#' @export
validate_logical <- function(logical){
    stopifnot("Please make sure all variables are logical" = all(is.logical(logical)))}

#' Validate min_count
#' 
#' This function validates a min_count.
#'
#' @param min_count A numeric.
#' @return A numeric.
#' @export
validate_min_count <- function(min_count){
    stopifnot("Please make sure that min_count is numeric and positive" = is.numeric(min_count) && min_count > 0)
    if(min_count < 1){
        message("min_count is less than 1, setting to 1")
        min_count <- 1}
    return(min_count)}

#' Validate pals
#' 
#' This function validates pals.
#'
#' @param dds A DESeq2 object.
#' @param group.by A character string of the group by.
#' @param pals A vector of colors.
#' @return A vector of colors.
validate_pals <- function(dds, group.by, pals){
    if(!is.null(pals)){
        stopifnot("Please make sure that length(pals) is equal or greater than the number of levels in group_by" = length(pals) >= length(levels(colData(dds)[[group.by]])))
        if(is.null(names(pals))){
            message("No names for pals, setting names to group_by levels")
            group.lv <- levels(colData(dds)[[group.by]])
            pals <- pals[1:length(group.lv)]
            names(pals) <- group.lv}
        stopifnot("Please make sure that the names(pals) contains all group_by levels" = all(names(pals) %in% levels(colData(dds)[[group.by]])))}
    return(pals)}

#' Validate shape
#' 
#' This function validates a shape.
#'
#' @param dds A DESeq2 object.
#' @param group.by A character string of the group by.
#' @param shape A character string of the shape.
#' @return A character string of the shape.
validate_shape <- function(dds, group.by, shape){
    if(!is.null(shape)){
        stopifnot("Please make sure that the shape is a column in colData" = shape %in% colnames(colData(dds)))}
    return(shape)}

#' Validate batch
#' 
#' This function validates a batch.
#'
#' @param dds A DESeq2 object.
#' @param batch A character string of the batch.
#' @return A character string of the batch.
validate_batch <- function(dds, batch){
    if(!is.null(batch)){
        stopifnot("Please make sure batch is a column in colData" = batch %in% colnames(colData(dds)))

        # check if batch is a factor
        if(!is.factor(colData(dds)[[batch]])){
            colData(dds)[[batch]] <- factor(colData(dds)[[batch]], levels = unique(colData(dds)[[batch]]))}

        # check if group_by is in the design
        if(class(dds) != "DESeqTransform"){
            des <- paste0(as.character(design(dds)), collapse = " ")
            if(!str_detect(des, batch)){
                message("batch is not in the design, adding to design")
                    design(dds) <- as.formula(paste0(des, " + ", batch))}
            }
        }
    return(dds)}

#' Validate order
#' 
#' This function validates an order.
#'
#' @param order A character string of the order.
#' @return A character string of the order.
validate_order <- function(order){
    stopifnot("Please make sure that the order is either pxfc, log2FoldChange, or padj" = order %in% c("pxfc", "log2FoldChange", "padj"))
    return(order)}

#' Validate msigdbr
#' 
#' This function validates msigdbr.
#'
#' @param msigdbr A character string of the msigdbr.
#' @return A character string of the msigdbr.
validate_msigdbr <- function(msigdbr){
    if(is.character(msigdbr)){
        if(str_detect(msigdbr, ".gmt$") & file.exists(msigdbr)){
            message("Reading msigdbr from specified .gmt file")
            msigdbr <- read_gmt(msigdbr) %>%
                pivot_longer(everything(), names_to = "gs_name", values_to = "gene_symbol") %>%
                filter(gene_symbol != "") %>%
                mutate(collection = gsub(".gmt", "", basename(msigdbr)))}}
    stopifnot("Please make sure that the msigdbr is a data frame and contains gs_name, gene_symbol, and collection columns" = is.data.frame(msigdbr) & all(c("gs_name", "gene_symbol", "collection") %in% colnames(msigdbr)))
    return(msigdbr)}

#' Validate gsea
#' 
#' This function validates gsea.
#'
#' @param gsea A data frame of gsea.
#' @return A data frame of gsea.
validate_gsea <- function(gsea){
    stopifnot("Please make sure that the gsea is a data frame and contains ID, NES, pvalue, qvalue, comparison, and collection columns" = is.data.frame(gsea) & all(c("ID", "NES", "pvalue", "qvalue", "comparison", "collection") %in% colnames(gsea)))
    for(c in c("NES", "pvalue", "qvalue")) {
        message(paste0("Converting ", c, " to numeric"))
        gsea[[c]] <- as.numeric(gsea[[c]])}
    return(gsea)}

#' Validate gsea for comparison
#' 
#' This function validates gsea for comparison.
#'
#' @param gsea A data frame of gsea.
#' @return A data frame of gsea.
validate_gsea_comparison <- function(gsea){
    gsea <- validate_gsea(gsea)
    stopifnot("Please make sure that the gsea contains only one comparison pair" = length(unique(gsea$comparison)) == 1)
    return(gsea)}
    
#' Validate gsea object
#' 
#' This function validates a gsea object.
#'  
#' @param gsea.obj A gsea object.
#' @return A gsea object.
validate_gsea_object <- function(gsea.obj){
    gsea.obj@result <- validate_gsea(gsea.obj@result)
    return(gsea.obj)}

#' Validate gene set
#' 
#' This function validates a gene set.
#'
#' @param gsea.obj A gsea object.
#' @return A gsea object.
validate_gene_set <- function(gsea, gene_set){
    stopifnot("Please make sure that the gene set is in the gsea ID column" = gene_set %in% gsea$ID)
    return(gene_set)}

#' Validate paths
#' 
#' This function validates paths.
#'
#' @param paths A character string of the paths.
#' @return A character string of the paths.
validate_paths <- function(paths){
    stopifnot("Please make sure that all paths exist" = all(dir.exists(paths)))}

