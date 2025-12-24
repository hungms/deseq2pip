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
run_gsea_pip <- function(res, org, msigdbr = import_msigdbr(org), save_data = T, save_dir = getwd()){
 
    # validations
    #res <- validate_res(res)
    #org <- validate_org(org)
    #msigdbr <- validate_msigdbr(msigdbr)
    #validate_paths(group_save_dir)

    # get comparisons
    comparisons <- unique(res$comparison)

    # perform GSEA for each comparison
    message("\n########################################################\nPerforming GSEA Analysis\n########################################################")
    gsea.list <- list()

    for(i in seq_along(comparisons)){

        res.comparison <- res %>% filter(comparison == comparisons[i])
        gsea.list[[i]] <- run_gsea(res = res.comparison, org = org, msigdbr = msigdbr, save_data = TRUE, save_dir = paste0(save_dir, "/", comparisons[i]))
        
        # generate GSEA plots for each collection
        for(j in unique(gsea.list[[i]]$collection)){
            gsea.collection <- gsea.list[[i]] %>% filter(collection == j)
            plot_gsea_barplot(gsea.collection, save_dir = paste0(save_dir, "/", comparisons[i]))}}

    return(gsea.list)}


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
    #res <- validate_res_comparison(res)
    #validate_org(org)
    #msigdbr <- validate_msigdbr(msigdbr)
    #validate_paths(save_dir)

    # get comparison
    comparison <- unique(res$comparison)
    message("<<", comparison, ">>")

    # split msigdbr by collection
    collectionsplit <- msigdbr$collection
    msigdbr.collection <- msigdbr %>%
        dplyr::select(-c("collection")) %>%
        split(., collectionsplit)

    duplicated_genes <- res$gene[duplicated(res$gene)]
    if(length(duplicated_genes) > 0){
        message(paste0("\t- removing ", length(duplicated_genes), " duplicated genes for GSEA: ", paste0(duplicated_genes, collapse = ", ")))
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
        message("\t- running GSEA for ", names(msigdbr.collection)[i], "...")
        gsea.obj <- quiet(GSEA(de.order, TERM2GENE = msigdbr.collection[[i]], pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 1000))
        collection <- names(msigdbr.collection)[i]
        
        # Skip if result is empty or has 0 rows
        if(is.null(gsea.obj@result) || nrow(gsea.obj@result) == 0){
            message(paste0("\t- skipping collection '", collection, "' - no GSEA results found"))
            next}
        
        gsea.obj@result$collection <- collection
        gsea.obj@result$comparison <- comparison

        if(save_data){
            message("\t- saving ", collection, " results to RDS and TSV...")
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
    #gsea <- validate_gsea_comparison(gsea)
    #validate_paths(save_dir)

    # get comparison and collection
    comparison <- unique(gsea$comparison)
    collection <- unique(gsea$collection)
    gsea$ID <- sub(paste0("^", collection, "_"), "", gsea$ID)
    message("\t- generating top ", n, " ", collection, " barplots...")

    # get selected pathways
    quiet({
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
            ggprism::theme_prism(border = T) + 
            theme(
                axis.line.y = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank()) +
            ylim(c(-yrange, yrange)) + 
            geom_hline(yintercept = 0, color = "grey40", size = 0.5) + 
            scale_x_discrete(expand = c(0.05, 0.05)) + 
            coord_flip()

            print(p)

        # save plot
        if (save_plot) {
            plot_name <- paste0("gsea_", collection, "_top", n)
            if (signif) {
                plot_name <- paste0(plot_name, "_signif")}
            save_plot(p, plot_name = paste0(plot_name, "_barplot.pdf"), save_dir = save_dir, w = 9, h = 7)}
        return(p)
    })
}
