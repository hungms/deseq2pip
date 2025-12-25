

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
#' @aliases validate_dds_atac
#' @export
validate_atac <- function(dds){
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

#' @rdname validate_atac
#' @export
validate_dds_atac <- validate_atac

#' Validate DESeq2 object for RNA-seq
#'
#' This function validates a DESeq2 object for RNA-seq.
#'
#' @param dds A DESeq2 object.
#' @return A DESeq2 object.
#' @export          
validate_var <- function(dds, var){

    # Handle vector of vars
    if(length(var) > 1){
        for(v in var){
            dds <- validate_var(dds, v)}
        return(dds)}
    
    # check if var is a factor
    if(!is.factor(colData(dds)[[var]])){
        colData(dds)[[var]] <- factor(colData(dds)[[var]], levels = unique(colData(dds)[[var]]))}

    dds <- validate_dds(dds)
    return(dds)}

#' Validate DESeq2 object for comparison
#'
#' This function validates a DESeq2 object for comparison.
#'
#' @param dds A DESeq2 object.
#' @return A DESeq2 object.
#' @export
validate_comparison <- function(dds, var, comparison){
    stopifnot("Please make sure that each group variable is present in var" = all(str_split(comparison, "_vs_") %>% unlist(.) %in% colData(dds)[[var]]))
    dds <- validate_var(dds, var)
    return(dds)}

#' Validate design formula
#'
#' This function validates a design formula and ensures all colData columns used in the design are factors with levels set.
#'
#' @param dds A DESeq2 object.
#' @param design A design formula (character string or formula object).
#' @return A DESeq2 object with validated and converted factors.
#' @export
validate_design <- function(dds, design){
    
    # Convert design to formula if it's a character string
    if(is.character(design)){
        design <- as.formula(design)}
    
    # Extract variable names from the design formula
    # all.vars() extracts all variable names from a formula
    vars_in_design <- all.vars(design)
    
    # Remove intercept term if present
    vars_in_design <- vars_in_design[vars_in_design != "."]
    
    # Check that all variables exist in colData
    missing_vars <- vars_in_design[!vars_in_design %in% colnames(colData(dds))]
    stopifnot("Please make sure that all variables in the design formula exist in colData(dds)" = length(missing_vars) == 0)
    
    # Ensure all variables are factors with levels set
    for(var in vars_in_design){
        if(!is.factor(colData(dds)[[var]])){
            message(paste0("Variable '", var, "' in design is not a factor, converting to factor"))
            colData(dds)[[var]] <- factor(colData(dds)[[var]], levels = unique(colData(dds)[[var]]))}
        else{
            # Check if factor has levels set
            if(length(levels(colData(dds)[[var]])) == 0){
                message(paste0("Variable '", var, "' is a factor but has no levels, setting levels"))
                colData(dds)[[var]] <- factor(colData(dds)[[var]], levels = unique(colData(dds)[[var]]))}
        }
    }
    
    return(dds)}

#' Validate comparisons
#'
#' This function validates that all comparisons listed are found in dds[[var]].
#' If comparisons is NULL, it will automatically generate all possible comparisons.
#'
#' @param dds A DESeq2 object.
#' @param var Column name in colData(dds) to use for grouping.
#' @param comparisons Vector of comparisons to validate. If NULL, generates all possible comparisons.
#' @return A validated vector of comparisons.
#' @export
validate_comparisons <- function(dds, var, comparisons){
    
    # Ensure var exists and is a factor
    stopifnot("Please make sure that var exists in colData(dds)" = var %in% colnames(colData(dds)))
    
    if(!is.factor(colData(dds)[[var]])){
        message("var is not a factor, converting to factor")
        colData(dds)[[var]] <- factor(colData(dds)[[var]], levels = unique(colData(dds)[[var]]))}
    
    # Generate comparisons if NULL
    if(is.null(comparisons)){
        comparisons <- generate_comparisons(levels(colData(dds)[[var]]))}
    
    # Get all valid levels from dds[[var]]
    valid_levels <- levels(colData(dds)[[var]])
    
    # Validate each comparison
    for(i in seq_along(comparisons)){
        comparison <- comparisons[i]
        
        # Split comparison by "_vs_"
        groups <- str_split(comparison, "_vs_") %>% unlist(.)
        
        # Handle group comparisons (with +)
        all_groups <- c()
        for(group in groups){
            if(str_detect(group, "\\+")){
                # Split by + and add individual groups
                sub_groups <- str_split(group, "\\+") %>% unlist(.) %>% str_trim(.)
                all_groups <- c(all_groups, sub_groups)}
            else{
                all_groups <- c(all_groups, group)}}
        
        # Check that all groups are in valid_levels
        missing_groups <- all_groups[!all_groups %in% valid_levels]
        if(length(missing_groups) > 0){
            stop(paste0("please make sure that all groups in comparison '", comparison, "' are present in dds[[var]]. Missing groups: ", paste0(missing_groups, collapse = ", ")))}
    }
    
    return(comparisons)}

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
        res[[c]] <- as.numeric(res[[c]])}
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
#' @param pals A vector of colors. If NULL, returns NULL.
#' @return A vector of colors or NULL.
#' @export
validate_pals <- function(dds, group.by, pals){
    if(is.null(pals)){
        return(pals)}
    
    stopifnot("Please make sure that length(pals) is equal or greater than the number of levels in var" = length(pals) >= length(levels(colData(dds)[[group.by]])))
    if(is.null(names(pals))){
        message("No names for pals, setting names to var levels")
        group.lv <- levels(colData(dds)[[group.by]])
        pals <- pals[1:length(group.lv)]
        names(pals) <- group.lv}
    stopifnot("Please make sure that the names(pals) contains all var levels" = all(names(pals) %in% levels(colData(dds)[[group.by]])))
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
#' @param batch A character string of the batch. If NULL, returns dds unchanged.
#' @return A DESeq2 object.
#' @export
validate_batch <- function(dds, batch){
    if(is.null(batch)){
        return(dds)}
    
    stopifnot("Please make sure batch is a column in colData" = batch %in% colnames(colData(dds)))

    # check if batch is a factor
    if(!is.factor(colData(dds)[[batch]])){
        colData(dds)[[batch]] <- factor(colData(dds)[[batch]], levels = unique(colData(dds)[[batch]]))}

    # check if var is in the design
    if(class(dds) != "DESeqTransform"){
        des <- paste0(as.character(design(dds)), collapse = " ")
        if(!str_detect(des, batch)){
            message("batch is not in the design, adding to design")
                design(dds) <- as.formula(paste0(des, " + ", batch))}
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
validate_gsea_result <- function(gsea){
    stopifnot("Please make sure that the gsea is a data frame and contains ID, NES, pvalue, qvalue, comparison, and collection columns" = is.data.frame(gsea) & all(c("ID", "NES", "pvalue", "qvalue", "comparison", "collection") %in% colnames(gsea)))
    for(c in c("NES", "pvalue", "qvalue")) {
        gsea[[c]] <- as.numeric(gsea[[c]])}
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
