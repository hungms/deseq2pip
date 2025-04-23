#' Write class file pipeline for EnrichmentMap  
#' 
#' This function writes a class file for EnrichmentMap visualization.
#' 
#' @param dds DESeq2 object
#' @param group_by Group by column name
#' @param cls_name Name of the class file. Default is "dds_class.cls"
#' @param save_dir Directory where the formatted files will be saved. Default is current working directory
write_cls <- function(dds, group_by, cls_name = "dds_class.cls", save_dir = getwd()){

    # validations
    dds <- validate_dds_group_by(dds, group_by)
    validate_paths(save_dir)

    # get group levels
    group.lv <- rev(levels(dds[[group_by]]))
    ngroup <- length(group.lv)

    # fill class matrix
    add_pad <- function(row, length) {
        c(row, rep("", length - length(row)))}
    ncol <- ifelse(ncol(dds) < 4, 4, ncol(dds))
    row1 <- add_pad(c(ncol(dds), ngroup, 1), ncol)
    row2 <- add_pad(c("#", group.lv), ncol)
    row3 <- as.character(dds[[group_by]])

    # write class file
    class <- matrix("", nrow = 3, ncol = ncol)
    class[1,] <- row1
    class[2,] <- row2
    class[3,] <- row3
    write.table(class, file = paste0(save_dir, "/", cls_name), row.names = F, col.names = F, quote = F)}

#' Write class file pipeline for EnrichmentMap
#' 
#' This function writes a class file for EnrichmentMap visualization.
#' 
#' @param dds DESeq2 object
#' @param group_by Group by column name
#' @param group_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param save_dir Directory where the formatted files will be saved. Default is current working directory  
#' @export
write_cls_pip <- function(dds, group_by, group_dir = group_dir, save_dir = getwd()){

    # validations
    dds <- validate_dds_group_by(dds, group_by)
    validate_paths(save_dir)

    # get comparison
    comparison <- list.files(group_dir)
    comparison <- comparison[grepl("_vs_", comparison)]
    group1 <- gsub("_vs_.*", "", comparison)
    group2 <- gsub(".*_vs_", "", comparison)
    one_to_all <- ifelse(length(unique(group2)) == 1, T, F)

    # add pad
    add_pad <- function(row, length) {
        c(row, rep("", length - length(row)))}

    if(one_to_all){
        message("Writing class files for ONE-TO-ALL comparisons...")
        for(i in seq_along(comparison)){
            cls_name <- paste0(comparison[i], "_class.cls")
            dds.onetoall <- dds
            colData(dds.onetoall)[[group_by]] <- ifelse(colData(dds.onetoall)[[group_by]] == group1[i], group1[i], group2[i])
            colData(dds.onetoall)[[group_by]] <- factor(colData(dds.onetoall)[[group_by]], levels = c(group2[i], group1[i]))
            write_cls(dds = dds.onetoall, group_by = group_by, cls_name = cls_name, save_dir = save_dir)}}
    else {
        message("Writing class files for PAIRWISE comparisons...")
        write_cls(dds = dds, group_by = group_by, cls_name = "dds_class.cls", save_dir = save_dir)}}


#' Format Enrichment Map Data
#'
#' This function processes and formats GSEA results for visualization in EnrichmentMap.
#'
#' @param dds DESeq2 object
#' @param org The organism to use, either "human" or "mouse".
#' @param collection Vector of gene set collections to process. Default includes HALLMARK, GOBP, KEGG, and REACTOME
#' @param group_dir Parent directory containing comparison subdirectories. Default is current working directory
#' @param save_dir Directory where the formatted files will be saved. Default is current working directory
#' @param save_dir_name Name of the subdirectory to save files in. Default is "enrichmentmap"
#' @return None. Creates formatted files for EnrichmentMap visualization
#' @examples
#' \dontrun{
#' # Format all default collections
#' enrichmentmap_pip(dds, group_by = "group", group_dir = "results")
#' 
#' # Format specific collections
#' enrichmentmap_pip(dds, group_by = "group", group_dir = "results", collection = c("HALLMARK", "KEGG"))
#' }
#' @export
enrichmentmap_pip <- function(
    dds,
    org,
    group_by,
    collection = c("HALLMARK", "GOBP", "KEGG", "REACTOME"),
    group_dir, 
    save_dir = group_dir,
    save_dir_name = "enrichmentmap"){

    dds <- validate_dds_group_by(dds, group_by)
    org <- validate_org(org)
    validate_paths(c(group_dir, save_dir))

    message("Formatting for enrichmentmap...")

    # create save directory
    enrich_save_dir <- paste0(save_dir, "/", save_dir_name)
    dir.create(enrich_save_dir, showWarnings = F, recursive = T)   

    #exprs
    counts <- assay(dds) %>% as.data.frame(.)
    write.table(counts, file = paste0(enrich_save_dir, "/dds_counts.txt"), row.names = T, col.names = T, quote = F)

    #gmt
    files <- list.files(system.file("extdata", package = "deseq2pip"), full.names = T)
    files <- files[which(str_detect(files, paste0(org, "_msigdbr.gmt")))]
    f <- files[length(files)]
    file.copy(from = f, to = paste0(enrich_save_dir, "/", org, "_msigdbr.gmt"), overwrite = TRUE)

    #class
    write_cls_pip(dds = dds, group_by = group_by, group_dir = group_dir, save_dir = enrich_save_dir)

    # for each comparison
    comparison <- list.files(group_dir, full.names = T)
    comparison <- comparison[grepl("_vs_", comparison)]

    for(i in seq_along(comparison)){
        
        #diffexp
        res <- read.table(paste0(comparison[i], "/diffexp_DESeq2.tsv"), sep = "\t", header = T)
        stopifnot(c("gene", "rank") %in% colnames(res))
        res.rank <- res %>% 
            select(gene, rank) %>%
            filter(rank != 0) %>%
            filter(!is.na(rank)) %>%
            filter(rank != "NA") %>%
            arrange(desc(rank))
        save_tsv(res.rank, tsv_name = paste0(basename(comparison[i]), "_diffexp_DESeq2_rank.rnk"), save_dir = enrich_save_dir)

        #gsea
        files <- list.files(comparison[i], full.names = T)
        collection.pattern <- paste0(collection, collapse = "|")
        stopifnot(any(str_detect(files, collection.pattern)))
        files <- files[str_detect(files, collection.pattern)]
        files <- files[which(str_detect(files, ".tsv$"))]
        
        enrichmentmap.list <- list()
        for(j in seq_along(files)){
            gsea.df <- read.table(files[j], sep = "\t", header = T)
            gsea.df <- gsea.df %>% 
                mutate(
                    phenotype = ifelse(NES > 0, "+1", "-1")) %>%
                select(ID, Description, pvalue, qvalue, phenotype)
            collection_name <- gsub("gsea_|\\.tsv", "", basename(files[j]))
            enrichmentmap.list[[j]] <- gsea.df}

        enrichmentmap.merged <- bind_rows(enrichmentmap.list)
        save_tsv(enrichmentmap.merged, tsv_name = paste0(basename(comparison[i]), "_enrichments.tsv"), save_dir = enrich_save_dir)}
    }