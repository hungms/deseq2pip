

#' Run Gene Set Enrichment Analysis
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using clusterProfiler
#' to identify enriched gene sets in the differential expression results. It can use
#' either pre-defined MSigDB gene sets or custom gene sets.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param org The organism to use, either "human" or "mouse". Default is "human"
#' @param group_by Column name in colData(dds) to use for grouping. Default is "Group1"
#' @param custom_msigdb Path to a custom gene set database in TSV format. Default is NULL
#' @param order Column name to use for ranking genes. Default is "rank"
#' @param save_data Logical. If TRUE, saves results to TSV and RDS files. Default is TRUE
#' @param save_dir Directory to save the results. Default is current working directory
#' @return A list of GSEA results for each gene set collection, with each element containing:
#'         - gene set name
#'         - normalized enrichment score (NES)
#'         - p-value and adjusted p-value
#'         - leading edge genes
#'         - collection name
#'         - comparison name
#' @examples
#' # Run GSEA with default MSigDB gene sets
#' gsea_results <- run_gsea(res)
#' 
#' # Run GSEA with custom gene sets
#' gsea_results <- run_gsea(res, custom_msigdb = "path/to/custom_sets.tsv")
#' 
#' # Run GSEA for mouse data
#' gsea_results <- run_gsea(res, org = "mouse")
#' @export
run_gsea <- function(res, org = "human", group_by = "Group1", custom_msigdb = NULL, order = "rank", save_data = T, save_dir = getwd()){
    if(length(custom_msigdb) > 0){
        stopifnot(str_detect(custom_msigdb, ".tsv$") & file.exists(custom_msigdb))
        msigdbr <- read.table(custom_msigdb, sep = "\t", header = T)}
    else{
        msigdbr <- import_msigdb(org = org)}

    stopifnot(c("gs_name", "gene_symbol", "collection") %in% colnames(msigdbr))
    stopifnot(c("padj", "gene", "log2FoldChange", "comparison") %in% colnames(res))
    stopifnot(order %in% c("rank", colnames(res)))
    
    comparison <- unique(res$comparison)
    stopifnot(length(comparison) == 1)

    res <- res %>% 
        mutate(rank = -log10(padj)*log2FoldChange)

    collectionsplit <- msigdbr$collection

    msigdbr.collection <- msigdbr %>%
        dplyr::select(-c("collection")) %>%
        split(., collectionsplit)

    de.order <- res %>%
        arrange(desc(.data[[order]])) %>%
        .$log2FoldChange

    names(de.order) <- res %>%
        arrange(desc(.data[[order]])) %>%
        .$gene
    
    de.order <- na.omit(de.order)
    de.order <- de.order[which(de.order != 0)]
    de.order = sort(de.order, decreasing = TRUE)

    gsea.list <- list()
    for(i in seq_along(msigdbr.collection)){

        gsea.obj <- GSEA(de.order, TERM2GENE = msigdbr.collection[[i]], pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 1000)
        collection <- names(msigdbr.collection)[i]
        comparison_name <- paste0(comparison, "_", collection)
        gsea.obj@result$collection <- collection
        gsea.obj@result$comparison <- comparison

        if(save_data){
            saveRDS(gsea.obj, file = paste0(save_dir, "/", comparison, "/gsea_", collection,".rds"))
            save_tsv(gsea.obj@result, tsv_name = paste0("gsea_", collection,".tsv"), save_dir = paste0(save_dir, "/", comparison, "/"))
            }

        gsea.list[[i]] <- gsea.obj
        names(gsea.list)[i] <- comparison_name}

    return(gsea.list)}



#' Run GSEA for All Comparisons
#'
#' This function performs Gene Set Enrichment Analysis for all comparisons in the
#' differential expression results.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @param ... Additional arguments passed to run_gsea()
#' @return A list of GSEA results for all comparisons and gene set collections
#' @examples
#' # Run GSEA for all comparisons
#' gsea_results <- run_gsea_list(res)
#' 
#' # Run with custom parameters
#' gsea_results <- run_gsea_list(res, org = "mouse")
#' @export
run_gsea_list <- function(res = NULL, save_dir = getwd(), ...){
    message("Running gene set enrichment analysis...")

    if(length(res) == 0){
        res <- read_diffexp_list(data_dir = save_dir)}

    stopifnot(is.data.frame(res))
    stopifnot(c("padj", "gene", "log2FoldChange", "comparison") %in% colnames(res))
    res.list <- split(res, res$comparison)

    gsea.list <- list()
    pb <- txtProgressBar(min = 0, max = length(res.list), style = 3)
    
    for(i in seq_along(res.list)){
        gsea.list <- c(gsea.list, run_gsea(res.list[[i]], save_dir = save_dir, ...))
        setTxtProgressBar(pb, i)
    }

    close(pb)
    return(gsea.list)
}



#' Read GSEA Results from RDS Files
#'
#' This function reads GSEA results from RDS files for multiple comparisons and collections.
#'
#' @param data_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param collection Vector of gene set collections to load. Default includes HALLMARK, GOBP, KEGG, REACTOME, BIOCARTA, and TFT
#' @return A list of GSEA objects for each comparison
#' @examples
#' # Read all default collections
#' gsea_results <- read_gsea_rds_list("results")
#' 
#' # Read specific collections
#' gsea_results <- read_gsea_rds_list("results", collection = c("HALLMARK", "KEGG"))
#' @export
read_gsea_rds_list <- function(
    data_dir = getwd(), 
    collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT")){
    
    comparison <- list.files(data_dir)
    files <- c()
    for(i in seq_along(comparison)){
        files <- c(files, paste0(comparison[i], "/", list.files(paste0(data_dir, "/", comparison[i]))))}
    files <- files[which(str_detect(files, ".rds"))]
    collection.pattern <- paste0(collection, collapse = "|")
    stopifnot(any(str_detect(files, collection.pattern)))
    files <- files[str_detect(files, collection.pattern)]

    gsea.list <- list()

    for(i in seq_along(files)){
        gsea.list[[i]] <- readRDS(paste0(data_dir, "/", files[i]))
        comparison <- gsub("/gsea", "", files[i])
        comparison <- gsub(".rds", "", comparison)
        names(gsea.list)[i] <- comparison}

    return(gsea.list)
}


#' Read GSEA Results from TSV Files
#'
#' This function reads GSEA results from TSV files for multiple comparisons and collections.
#'
#' @param data_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param collection Vector of gene set collections to load. Default includes HALLMARK, GOBP, KEGG, REACTOME, BIOCARTA, and TFT
#' @param merge Logical. If TRUE, returns a single merged data frame. If FALSE, returns a list of data frames
#' @return Either a merged data frame or a list of data frames containing GSEA results
#' @examples
#' # Read and merge all results
#' merged_gsea <- read_gsea_tsv_list("results", merge = TRUE)
#' 
#' # Read specific collections as a list
#' gsea_list <- read_gsea_tsv_list("results", collection = c("HALLMARK", "KEGG"), merge = FALSE)
#' @export
read_gsea_tsv_list <- function(
    data_dir = getwd(), 
    collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT"),
    merge = T){

    comparison <- list.files(data_dir)
    files <- c()
    for(i in seq_along(comparison)){
        files <- c(files, paste0(comparison[i], "/", list.files(paste0(data_dir, "/", comparison[i]))))}
    files <- files[which(str_detect(files, ".tsv"))]
    files <- files[which(!str_detect(files, "enrichmentmap"))]
    collection.pattern <- paste0(collection, collapse = "|")
    stopifnot(any(str_detect(files, collection.pattern)))
    files <- files[str_detect(files, collection.pattern)]
        
    gsea.df <- list()

    for(i in seq_along(files)){
        gsea.df[[i]] <- read.table(paste0(data_dir, "/", files[i]), sep = "\t", header = T)}

    if(merge){
        gsea.df <- bind_rows(gsea.df)}

    return(gsea.df)
}

#' Format Enrichment Map Data
#'
#' This function processes and formats GSEA results for visualization in EnrichmentMap.
#'
#' @param data_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param collection Vector of gene set collections to process. Default includes HALLMARK, GOBP, KEGG, and REACTOME
#' @param save_dir Directory where the formatted files will be saved. Default is current working directory
#' @return None. Creates formatted files for EnrichmentMap visualization
#' @examples
#' # Format all default collections
#' format_enrichmentmap("results")
#' 
#' # Format specific collections
#' format_enrichmentmap("results", collection = c("HALLMARK", "KEGG"))
#' @export
format_enrichmentmap <- function(
    data_dir = getwd(), 
    collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME"),
    save_dir = data_dir){
    message("Formatting for enrichmentmap...")
    comparison <- list.files(data_dir)

    for(i in seq_along(comparison)){
        files <- list.files(paste0(data_dir, "/", comparison[i]))

        #diffexp
        res <- read.table(paste0(data_dir, "/", comparison[i], "/diffexp_deseq2_wald.tsv"), sep = "\t", header = T)
        res.rank <- res %>% 
            select(gene, rank)
        save_tsv(res.rank, tsv_name = paste0("diffexp_deseq2_wald_rank.rnk"), save_dir = paste0(save_dir, "/", comparison[i], "/enrichmentmap/"))

        #gsea
        collection.pattern <- paste0(collection, collapse = "|")
        stopifnot(any(str_detect(files, collection.pattern)))
        files <- files[str_detect(files, collection.pattern)]
        files  <- files[which(str_detect(files, ".tsv$"))]
        
        enrichmentmap.list <- list()
        for(j in seq_along(files)){
            gsea.df <- read.table(paste0(data_dir, "/", comparison[i], "/", files[j]), sep = "\t", header = T)
            gsea.df <- gsea.df %>% 
                mutate(
                    phenotype = ifelse(NES > 0, "+1", "-1")) %>%
                select(ID, Description, pvalue, qvalue, phenotype)
            collection_name <- gsub("gsea_|\\.tsv", "", files[j])
            save_tsv(gsea.df, tsv_name = paste0("gsea_", collection_name, "_enrichmentmap.tsv"), save_dir = paste0(save_dir, "/", comparison[i], "/enrichmentmap/"))
            enrichmentmap.list[[j]] <- gsea.df}

        enrichmentmap.merged <- bind_rows(enrichmentmap.list)
        save_tsv(enrichmentmap.merged, tsv_name = "gsea_enrichmentmap_merged.tsv", save_dir = paste0(save_dir, "/", comparison[i], "/enrichmentmap/"))
        }
        }


#' Generate GSEA Barplots for All Comparisons
#'
#' This function creates barplots to visualize GSEA results for all comparisons
#' and gene set collections.
#'
#' @param gsea.df GSEA result data frame from run_gsea()
#' @param ... Additional arguments passed to plot_gsea_barplot()
#' @return A list of GSEA barplots for each comparison and collection
#' @examples
#' # Generate GSEA barplots for all results
#' plots <- plot_gsea_barplot_list(gsea_results)
#' 
#' # Generate with custom parameters
#' plots <- plot_gsea_barplot_list(gsea_results, n = 15, signif = TRUE)
#' @export
plot_gsea_barplot_list <- function(gsea.df, ...){
    message("Generating GSEA barplots...")
    stopifnot(is.data.frame(gsea.df))
    gsea.df.list <- split(gsea.df, interaction(gsea.df$comparison, gsea.df$collection))

    gsea.barplot.list <- list()
    pb <- txtProgressBar(min = 0, max = length(gsea.df.list), style = 3)
    
    for(i in seq_along(gsea.df.list)){
        gsea.barplot <- plot_gsea_barplot(gsea.df.list[[i]], ...)
        gsea.barplot.list <- c(gsea.barplot.list, gsea.barplot)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    return(gsea.barplot.list)
}