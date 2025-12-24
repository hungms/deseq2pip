#' Summarize Gene Expression in DESeq2 Object
#'
#' This function aggregates expression values from multiple isoforms/peaks of the same gene
#' in a DESeq2 object into a single gene-level expression value. It creates a new
#' DESeq2 object with gene-level expression data.
#'
#' @param dds DESeq2 object containing isoform/peak-level expression data
#' @param gene_sym_col Column name in rowData(dds) containing gene symbols. Default is "Gene.Name"
#' @param ... Additional arguments passed to summarize_genes()
#' @return A new DESeq2 object containing gene-level expression data
#' @examples
#' \dontrun{
#' # Sum gene expression
#' gene_dds <- summarize_genes_dds(dds)
#' 
#' # Average gene expression
#' gene_dds <- summarize_genes_dds(dds, normalized = TRUE)
#' }
#' @export
summarize_genes_dds <- function(dds, gene_sym_col = "gene", ...){
    df <- assay(dds) %>% as.data.frame(.)
    gene_sym_vec <- rowData(dds)[[gene_sym_col]]
    df <- strpip::summarize_genes(input = df, gene_sym_vec = gene_sym_vec, ...)

    metadata <- as.data.frame(colData(dds))[colnames(df),]
    dds.new <- DESeqDataSetFromMatrix(
        countData = df,
        colData = metadata,
        design = design(dds))
    rowData(dds.new)$gene <- rownames(df)
    return(dds.new)
}

#' Import nfcore/rnaseq DESeq2 Object
#'
#' This function imports a DESeq2 object from nfcore RNA-seq pipeline output and adds gene names.
#'
#' @param rdata Path to the RData file containing the DESeq2 object
#' @param tx2gene Path to the tx2gene mapping file
#' @return A DESeq2 object with gene names added to rowData
#' @examples
#' \dontrun{
#' # Import DESeq2 object and add gene names
#' dds <- import_nfcore_rna("path/to/dds.RData", "path/to/tx2gene.tsv")
#' }
#' @export
import_nfcore_rna <- function(rdata, tx2gene){
    # Load DESeq2 object
    load(rdata)
    dds <- get("dds")

    # Read gene symbol mapping and ensure it's properly formatted
    genes <- read.table(tx2gene, sep = "\t", header = F) %>%
        dplyr::mutate(V3 = as.character(V3)) %>%  # Explicitly convert gene symbols to character
        dplyr::distinct(V2, V3) %>%
        tibble::column_to_rownames("V2")
    
    # Ensure rownames of dds are character vectors
    rownames(dds) <- as.character(rownames(dds))
    rowData(dds)$gene <-  as.character(genes[rownames(dds),])

    # Summarize genes
    dds.summarized <- summarize_genes_dds(dds, gene_sym_col = "gene", normalize = F)
    return(dds.summarized)
}

#' Import nfcore ATAC-seq DESeq2 Object
#'
#' This function imports a DESeq2 object from nfcore ATAC-seq pipeline output and adds gene names.
#'
#' @param rdata Path to the RData file containing the DESeq2 object
#' @param tx2gene Path to the tx2gene mapping file
#' @param dist.to.TSS Distance to TSS to consider as TSS. Default is 2000
#' @return A DESeq2 object with gene names added to rowData
#' @examples
#' \dontrun{
#' # Import DESeq2 object and add gene names
#' dds <- import_nfcore_atac("path/to/dds.RData", "path/to/tx2gene.tsv")
#' }
#' @export
import_nfcore_atac <- function (rdata, annotatePeaks, dist.to.TSS = 2000) 
{
    load(rdata)
    dds <- get("dds")
    rowData(dds)$peaks <- rownames(dds)
    annotations <- read.table(annotatePeaks, sep = "\t", header = T, row.names = 1)
    annotations <- annotations[rownames(assay(dds)), ]
    rowData(dds) <- cbind(rowData(dds), annotations)
    rowData(dds)$gene <- rowData(dds)$Gene.Name
    rowData(dds)$Gene.Name <- NULL
    rowData(dds)$TSS <- ifelse(rowData(dds)$gene != "" & rowData(dds)$Distance.to.TSS > -dist.to.TSS & rowData(dds)$Distance.to.TSS < dist.to.TSS, T, F)
    return(dds)
}


#' Read Differential Expression Results
#'
#' This function reads differential expression results from multiple comparison directories
#' and optionally merges them into a single data frame.
#'
#' @param group_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param merge Logical. If TRUE, returns a single merged data frame. If FALSE, returns a list of data frames
#' @return Either a merged data frame or a list of data frames containing differential expression results
#' @examples
#' \dontrun{
#' # Read and merge all results
#' merged_results <- read_diffexp(group_dir = "results", merge = TRUE)
#' 
#' # Read results as a list
#' results_list <- read_diffexp(group_dir = "results", merge = FALSE)
#' }
#' @export
read_diffexp <- function(group_dir = getwd(), merge = T){

    # validations
    validate_logical(merge)
    validate_paths(group_dir)

    # set comparison
    comparison_vec <- list.files(group_dir)
    paths <- paste0(group_dir, "/", comparison_vec) # change to full path
    validate_paths(paths)

    # get diffexp files
    files <- c()
    for(i in seq_along(paths)){
        files <- c(files, list.files(paths[i], full.names = T))}
    files <- files[which(basename(files) == "diffexp_DESeq2.tsv")]

    # read diffexp files
    res.list <- list()
    for(i in seq_along(files)){

        res.list[[i]] <- read.table(files[i], sep = "\t", header = T)
        res.list[[i]] <- validate_res(res.list[[i]])

        # add comparison column if it does not exist
        if(!"comparison" %in% colnames(res.list[[i]])){
            comparison <- basename(gsub("diffexp_DESeq2.tsv", "", files[i]))
            res.list[[i]]$comparison <- comparison
            names(res.list)[i] <- comparison}
        else{
            names(res.list)[i] <- unique(res.list[[i]]$comparison)}}

    # merge diffexp files
    if(merge){
        res.list <- bind_rows(res.list)}

    return(res.list)
}


#' Read Differential Expression Results
#'
#' This function reads differential expression results from multiple comparison directories
#' and optionally merges them into a single data frame.
#'
#' @param group_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param merge Logical. If TRUE, returns a single merged data frame. If FALSE, returns a list of data frames
#' @return Either a merged data frame or a list of data frames containing differential expression results
#' @examples
#' \dontrun{
#' # Read and merge all results
#' merged_results <- read_diffexp(group_dir = "results", merge = TRUE)
#' 
#' # Read results as a list
#' results_list <- read_diffexp(group_dir = "results", merge = FALSE)
#' }
#' @export
read_diffexp <- function(group_dir = getwd(), merge = T){

    # validations
    validate_logical(merge)
    validate_paths(group_dir)

    # set comparison
    comparison_vec <- list.files(group_dir)
    paths <- paste0(group_dir, "/", comparison_vec) # change to full path
    validate_paths(paths)

    # get diffexp files
    files <- c()
    for(i in seq_along(paths)){
        files <- c(files, list.files(paths[i], full.names = T))}
    files <- files[which(basename(files) == "diffexp_DESeq2.tsv")]

    # read diffexp files
    res.list <- list()
    for(i in seq_along(files)){

        res.list[[i]] <- read.table(files[i], sep = "\t", header = T)
        res.list[[i]] <- validate_res(res.list[[i]])

        # add comparison column if it does not exist
        if(!"comparison" %in% colnames(res.list[[i]])){
            comparison <- basename(gsub("diffexp_DESeq2.tsv", "", files[i]))
            res.list[[i]]$comparison <- comparison
            names(res.list)[i] <- comparison}
        else{
            names(res.list)[i] <- unique(res.list[[i]]$comparison)}}

    # merge diffexp files
    if(merge){
        res.list <- bind_rows(res.list)}

    return(res.list)
}



#' Read GSEA Results from TSV Files
#'
#' This function reads GSEA results from TSV files for multiple comparisons and collections.
#'
#' @param group_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param collection Vector of gene set collections to load. Default includes HALLMARK, GOBP, KEGG, REACTOME, BIOCARTA, and TFT
#' @param merge Logical. If TRUE, returns a single merged data frame. If FALSE, returns a list of data frames
#' @return Either a merged data frame or a list of data frames containing GSEA results
#' @examples
#' \dontrun{
#' # Read and merge all results
#' merged_gsea <- read_gsea_tsv("results", merge = TRUE)
#' 
#' # Read specific collections as a list
#' gsea_list <- read_gsea_tsv("results", collection = c("HALLMARK", "KEGG"), merge = FALSE)
#' }
#' @export
read_gsea_tsv <- function(
    group_dir = getwd(), 
    collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT"),
    merge = T){

    # validations
    validate_paths(group_dir)

    # get comparison
    comparison <- list.files(group_dir)
    files <- c()
    for(i in seq_along(comparison)){
        files <- c(files, paste0(comparison[i], "/", list.files(paste0(group_dir, "/", comparison[i]))))}
    files <- files[which(str_detect(files, ".tsv"))]
    files <- files[which(!str_detect(files, "enrichmentmap"))]
    collection.pattern <- paste0(collection, collapse = "|")
    stopifnot(any(str_detect(files, collection.pattern)))
    files <- files[str_detect(files, collection.pattern)]
        
    gsea <- list()

    for(i in seq_along(files)){
        gsea[[i]] <- read.table(paste0(group_dir, "/", files[i]), sep = "\t", header = T)}

    if(merge){
        gsea <- bind_rows(gsea)}

    return(gsea)
}


#' Read GSEA Results from RDS Files
#'
#' This function reads GSEA results from RDS files for multiple comparisons and collections.
#'
#' @param group_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param collection Vector of gene set collections to load. Default includes HALLMARK, GOBP, KEGG, REACTOME, BIOCARTA, and TFT
#' @return A list of GSEA objects for each comparison
#' @examples
#' \dontrun{
#' # Read all default collections
#' gsea_results <- read_gsea_rds("results")
#' 
#' # Read specific collections
#' gsea_results <- read_gsea_rds("results", collection = c("HALLMARK", "KEGG"))
#' }
#' @export
read_gsea_rds <- function(
    group_dir = getwd(), 
    collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT")){

    # validations
    validate_paths(group_dir)

    # get comparison
    comparison <- list.files(group_dir)
    files <- c()
    for(i in seq_along(comparison)){
        files <- c(files, paste0(comparison[i], "/", list.files(paste0(group_dir, "/", comparison[i]))))}
    files <- files[which(str_detect(files, ".rds"))]
    collection.pattern <- paste0(collection, collapse = "|")
    stopifnot(any(str_detect(files, collection.pattern)))
    files <- files[str_detect(files, collection.pattern)]

    gsea.list <- list()

    for(i in seq_along(files)){
        gsea.list[[i]] <- readRDS(paste0(group_dir, "/", files[i]))
        comparison <- gsub("/gsea", "", files[i])
        comparison <- gsub(".rds", "", comparison)
        names(gsea.list)[i] <- comparison}

    return(gsea.list)
}