#' Import MSigDB Gene Sets
#'
#' This function imports pre-defined MSigDB gene sets for human or mouse organisms.
#'
#' @param org The organism to use, either "human" or "mouse".
#' @return A data frame containing MSigDB gene sets
#' @examples
#' \dontrun{
#' # Import human MSigDB gene sets
#' msigdbr_human <- import_msigdbr("human")
#' 
#' # Import mouse MSigDB gene sets
#' msigdbr_mouse <- import_msigdbr("mouse")
#' }
#' @export
import_msigdbr <- function(org) {
    files <- list.files(system.file("extdata", package = "deseq2pip"), full.names = T)
    files <- files[which(str_detect(files, paste0(org, "_msigdbr.tsv")))]
    f <- files[length(files)]
    msigdbr <- read.table(gzfile(f), header = TRUE, sep = "\t") %>% as.data.frame(.)
    return(msigdbr)}

#' Run Gene Set Enrichment Analysis
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using clusterProfiler
#' to identify enriched gene sets in the differential expression results. It can use
#' either pre-defined MSigDB gene sets or custom gene sets.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param org The organism to use, either "human" or "mouse".
#' @param custom_msigdb Dataframe or path to a custom gene set database in GMT format. Default is NULL
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
#' \dontrun{
#' # Run GSEA with default MSigDB gene sets
#' gsea_results <- run_gsea(res)
#' 
#' # Run GSEA with custom gene sets
#' gsea_results <- run_gsea(res, custom_msigdb = "path/to/custom_sets.tsv")
#' 
#' # Run GSEA for mouse data
#' gsea_results <- run_gsea(res, org = "mouse")
#' }
#' @export
run_gsea <- function(res, org, msigdbr = import_msigdbr(org), save_data = T, save_dir = getwd()){

    # validations
    res <- validate_res_comparison(res)
    validate_org(org)
    msigdbr <- validate_msigdbr(msigdbr)
    validate_paths(save_dir)

    # split msigdbr by collection
    collectionsplit <- msigdbr$collection
    msigdbr.collection <- msigdbr %>%
        dplyr::select(-c("collection")) %>%
        split(., collectionsplit)

    # select comparison
    comparison <- unique(res$comparison)
    duplicated_genes <- res$gene[duplicated(res$gene)]
    if(length(duplicated_genes) > 0){
        message(paste0("Removing ", length(duplicated_genes), " duplicated genes for GSEA: ", paste0(duplicated_genes, collapse = ", ")))
        res <- res %>%
            filter(!gene %in% duplicated_genes)}

    # get gene rankings
    de.order <- res %>%
        arrange(desc(.data[["rank"]])) %>%
        .$rank

    # get gene names
    names(de.order) <- res %>%
        arrange(desc(.data[["rank"]])) %>%
        .$gene
    
    # remove 0s
    de.order <- na.omit(de.order)
    de.order <- de.order[which(de.order != 0)]
    de.order = sort(de.order, decreasing = TRUE)

    # run gsea for each collecition
    gsea.list <- list()
    for(i in seq_along(msigdbr.collection)){

        gsea.obj <- GSEA(de.order, TERM2GENE = msigdbr.collection[[i]], pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 1000)
        collection <- names(msigdbr.collection)[i]
        gsea.obj@result$collection <- collection
        gsea.obj@result$comparison <- comparison

        if(save_data){
            saveRDS(gsea.obj, file = paste0(save_dir, "/gsea_", collection,".rds"))
            save_tsv(gsea.obj@result, tsv_name = paste0("gsea_", collection,".tsv"), save_dir = save_dir)}

        gsea.list[[i]] <- gsea.obj@result
        names(gsea.list)[i] <- collection}
    
    gsea <- bind_rows(gsea.list)
    return(gsea)}


#' Plot GSEA Results
#'
#' This function creates a barplot to visualize GSEA results, showing the most significant
#' gene sets in each direction (up and down-regulated).
#'
#' @param gsea GSEA result data frame from run_gsea()
#' @param n Number of top gene sets to show in each direction. Default is 10
#' @param signif Logical. If TRUE, only shows gene sets with q-value < 0.05. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is current working directory
#' @return A ggplot object showing the GSEA barplot
#' @examples
#' \dontrun{
#' # Basic GSEA plot
#' p <- plot_gsea_barplot(gsea_results)
#' 
#' # GSEA plot with more gene sets and only significant results
#' p <- plot_gsea_barplot(gsea_results, n = 15, signif = TRUE)
#' }
#' @export
plot_gsea_barplot <- function(gsea, n = 10, signif = F, save_plot = T, save_dir = getwd()){

    # validations
    gsea <- validate_gsea_comparison(gsea)
    validate_paths(save_dir)

    # get comparison and collection
    comparison <- unique(gsea$comparison)
    collection <- unique(gsea$collection)
    gsea$ID <- sub(paste0("^", collection, "_"), "", gsea$ID)

    # get selected pathways
    selected_pathways <- gsea %>% 
        arrange(desc(NES^2)) %>% 
        mutate(direction = ifelse(NES > 0, "Up", "Down"), direction = factor(direction, c("Up", "Down"))) %>% 
        group_by(direction) %>% 
        slice_min(n = n, order_by = qvalue, with_ties = F) %>% 
        .$ID

    # get significant pathways
    if (signif){
        selected_pathways <- gsea %>% 
            filter(qvalue < 0.05) %>% 
            arrange(desc(NES^2)) %>% 
            mutate(direction = ifelse(NES > 0, "Up", "Down"), direction = factor(direction, c("Up", "Down"))) %>% 
            group_by(direction) %>% 
            slice_min(n = n, order_by = qvalue, with_ties = F) %>% 
            .$ID}

    # get selected gsea
    selected.gsea <- gsea %>% 
        filter(ID %in% selected_pathways) %>% 
        select(ID, NES, pvalue, qvalue)

    # get empty row
    empty_row <- data.frame(ID = "", NES = 0, pvalue = 1, qvalue = 1)
    yrange <- max(abs(selected.gsea$NES)) * 1.1

    # plot
    p <- selected.gsea %>% 
        rbind(., empty_row) %>% 
        mutate(label = case_when(
            qvalue < 0.001 ~ "***", 
            qvalue < 0.01 ~ "**", 
            qvalue < 0.05 ~ "*", 
            .default = "")) %>% 
        mutate(size = ifelse(nchar(ID) >= 50, 1, ifelse(nchar(ID) >= 35, 2, 3)), 
            ID = str_replace_all(ID, ".{50}", "\\0\n")) %>% 
        ggplot(aes(x = fct_reorder(ID, NES), y = NES, fill = NES)) + 
        geom_col(aes(stroke = label), width = 0.75, col = "black") + 
        geom_text(aes(y = ifelse(NES > 0, -0.1, 0.1), label = fct_reorder(ID, NES), hjust = ifelse(NES > 0, 1, 0), size = size), fontface = "bold") + 
        scale_size(range = c(2, 3)) + 
        geom_text(aes(label = label, y = NES + 0.2 * sign(NES)), 
            position = position_dodge(width = 0.75), vjust = 0.75, 
            size = 5) + 
        xlab(NULL) + 
        scale_fill_distiller(palette = "RdBu") + 
        guides(
            fill = guide_colorbar(
                title = "NES", 
                title.position = "top", 
                direction = "vertical", 
                frame.colour = "black", 
                ticks.colour = "black", 
                order = 1), 
            size = guide_none()) + 
        ggtitle(comparison, subtitle = paste0(collection, " COLLECTION")) + 
        theme_border() + 
        theme_text() + 
        theme(
            panel.border = element_rect(fill = NA, color = "black", size = 0.7), 
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
            plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5), 
            axis.line.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank()) + 
        theme_gridlines() + 
        ylim(c(-yrange, yrange)) + 
        geom_hline(yintercept = 0, color = "grey40", size = 0.5) + 
        scale_x_discrete(expand = c(0.05, 0.05)) + 
        coord_flip()

    # save plot
    if (save_plot) {
        plot_name <- paste0("gsea_", collection, "_top", n)
        if (signif) {
            plot_name <- paste0(plot_name, "_signif")}
        save_plot(p, plot_name = paste0(plot_name, "_barplot.pdf"), save_dir = save_dir, w = 9, h = 7)}
    print(p)
    return(p)
}

#' Run GSEA wrapper
#'
#' This function runs GSEA for all comparisons and collections.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param org The organism to use, either "human" or "mouse".
#' @param msigdbr Dataframe or path to a custom gene set database in GMT format. Default is NULL
#' @param save_data Logical. If TRUE, saves results to RDS and TSV files. Default is TRUE   
#' @param group_save_dir Directory to save the results.
#' @return A list of GSEA results for each gene set collection
#' @examples
#' \dontrun{
#' # Run GSEA wrapper
#' gsea_results <- run_gsea_wrapper(res)
#' }
#' @export
run_gsea_wrapper <- function(res, org, msigdbr = import_msigdbr(org), save_data = T, group_save_dir){

    # validations
    res <- validate_res_comparison(res)
    org <- validate_org(org)
    msigdbr <- validate_msigdbr(msigdbr)
    validate_paths(group_save_dir)

    # get comparison
    comparison <- unique(res$comparison)
    comparison_save_dir <- paste0(group_save_dir, "/", comparison)
    validate_paths(comparison_save_dir)

    # Run GSEA
    gsea <- run_gsea(res = res, org = org, msigdbr = msigdbr, save_dir = comparison_save_dir)
    
    # Generate GSEA plots for each collection
    for(i in unique(gsea$collection)){
        gsea.collection <- gsea %>% filter(collection == i)
        plot_gsea_barplot(gsea.collection, save_dir = comparison_save_dir)}
    
    return(gsea)}


#' Run GSEA pipeline
#' 
#' This function runs the GSEA pipeline for all comparisons and collections.
#' 
#' @param res Differential expression result data frame from run_diffexp()
#' @param org The organism to use, either "human" or "mouse".
#' @param msigdbr Dataframe or path to a custom gene set database in GMT format. Default is NULL
#' @param save_data Logical. If TRUE, saves results to RDS and TSV files. Default is TRUE
#' @param group_save_dir Directory to save the results.
#' @return A list of GSEA results for each gene set collection
#' @export 
run_gsea_pip <- function(res, org, msigdbr = import_msigdbr(org), save_data = T, group_save_dir){
 
    # validations
    res <- validate_res(res)
    org <- validate_org(org)
    msigdbr <- validate_msigdbr(msigdbr)
    validate_paths(group_save_dir)

    # generate a vector of comparison
    comparison_vec <- unique(res$comparison)
    validate_paths(group_save_dir)

    # run gsea for each comparison
    gsea.list <- list()
    for(i in seq_along(comparison_vec)){
        res.comparison <- res %>% filter(comparison == comparison_vec[i])
        gsea.list[[i]] <- run_gsea_wrapper(res = res.comparison, org = org, msigdbr = msigdbr, group_save_dir = group_save_dir)}

    # bind rows
    gsea <- bind_rows(gsea.list)

    return(gsea)}


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