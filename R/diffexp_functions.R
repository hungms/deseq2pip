#' Generate Comparison Names
#' 
#' This function generates comparison names for a given vector of group levels.
#' It takes a vector of group levels and returns a vector of comparison names.
#' 
#' @param group_vec A vector of group levels
#' @return A vector of comparison names
#' @export
generate_comparisons <- function(group_vec){

    # generate comparisons
    combos <- combn(group_vec, 2)
    comparisons_vec <- apply(combos, 2, function(pair) {
        if (pair[2] < pair[1]) paste0(c(as.character(pair[2]), as.character(pair[1])), collapse = "_vs_")
        else paste0(c(as.character(pair)), collapse = "_vs_")})    
    # return comparison names
    return(comparisons_vec)
}

#' Run Differential Expression Analysis
#'
#' This function performs differential expression analysis using DESeq2 for a comparison
#' between two groups. It calculates fold changes, p-values, and adjusted p-values for
#' each gene, and ranks genes based on their statistical and biological significance.
#'
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse".
#' @param group_by Column name in colData(dds) to use for grouping.
#' @param comparison Comparison to run.
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param save_data Logical. If TRUE, saves results to TSV file. Default is TRUE
#' @param save_dir Directory to save the results. Default is current working directory
#' @param ... Additional arguments passed to DESeq2::DESeq()
#' @return A data frame containing differential expression results with columns:
#'         - gene: gene identifier
#'         - comparison: name of the comparison
#'         - log2FoldChange: log2 fold change between groups
#'         - padj: adjusted p-value
#'         - pxfc: combined score (-log10(padj) * log2FoldChange)
#' @examples
#' \dontrun{
#' # Basic differential expression analysis
#' res <- run_diffexp(dds)
#' 
#' # Differential expression analysis with custom parameters
#' res <- run_diffexp(dds, group_by = "Treatment", save_data = TRUE)
#' }
#' @export
run_diffexp <- function(dds, org, group_by, comparison, order = "pxfc", save_data = T, save_dir = getwd(), ...){

    # validations
    dds <- validate_dds_comparison(dds, group_by, comparison)
    org <- validate_org(org)
    order <- validate_order(order)
    validate_paths(save_dir)
    org.to <- ifelse(org == "human", "mouse", "human")

    # subset dds
    comparison_groups <- str_split(comparison, "_vs_") %>% unlist(.)
    dds <- dds[,which(colData(dds)[[group_by]] %in% comparison_groups)]
    dds[[group_by]] <- droplevels(dds[[group_by]])

    # run deseq2
    dds <- DESeq(dds, ...)
    # Explicitly specify contrast: numerator vs reference
    # For "A_vs_B": A (numerator, comparison_groups[1]) vs B (reference, comparison_groups[2])
    # This gives positive log2FC when A > B, which matches the comparison name
    contrast_vec <- c(group_by, comparison_groups[1], comparison_groups[2])
    res <- results(dds, contrast = contrast_vec)

    # shrink lfc by ashr (contrast is already specified in res object)
    res <- lfcShrink(dds, res = res, type = "ashr") %>% 
        as.data.frame(.)
    
    # add existing row metadata to results
    row.meta <- rowData(dds) %>% as.data.frame(.)
    row.meta[which(colnames(row.meta) %in% colnames(res))] <- NULL
    if(ncol(row.meta) > 0){
        message("Combining row metadata to results...")
        if(all(rownames(res) == rownames(assay(dds)))){
            res <- res[rownames(assay(dds)),]
            res <- cbind(res, row.meta)}}

    # add pxfc and rank columns
    res <- res %>%
        mutate(
            comparison = comparison,
            pxfc = -log10(padj)*log2FoldChange,
            pxfc = ifelse(is.infinite(pxfc), 0, pxfc),
            rank = !!sym(order))

    # annotate and filter results
    res <- run_annotation(res, org.from = org, org.to = org.to, gene_column = "gene") %>%
        filter(!is.na(padj)) %>%
        filter(padj != "NA") %>%
        arrange(desc(rank))

    # save data
    if(save_data){
        save_tsv(res, tsv_name = "diffexp_DESeq2.tsv", save_dir = save_dir)}

    return(res)}

#' Generate MA Plot
#'
#' This function creates a MA plot to visualize differential expression results.
#' It shows the relationship between mean expression and log2 fold change for each gene.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param save_plot Logical. If TRUE, saves the plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is current working directory
#' @return A ggplot object showing the MA plot
#' @examples
#' \dontrun{
#' # Basic MA plot
#' p <- plot_ma(res)
#' }
#' @export
plot_ma <- function(res, fc.thresh = 0.5, save_plot = TRUE, save_dir = getwd()) {

    # validations
    res <- validate_res_comparison(res)
    validate_paths(save_dir)
    comparison <- unique(res$comparison)


    # set plot parameters
    res <- res %>% 
        mutate(baseMean = log2(baseMean + 1)) %>% 
        filter(baseMean >= 3) %>% 
        mutate(direction = case_when(
            log2FoldChange > fc.thresh & padj < 0.05 ~ "Up", 
            log2FoldChange < -fc.thresh & padj < 0.05 ~ "Down", 
            .default = "Non-DE"), 
            direction = factor(direction, c("Up", "Down", "Non-DE")), 
            size = ifelse(direction == "Non-DE", 0.5, 1.5))

    # set filtering parameters
    ngene <- nrow(res)
    res.label <- res %>% filter(direction != "Non-DE")
    mean.upper <- quantile(res.label$baseMean, probs = 0.999)
    fc.upper <- quantile(abs(res.label$log2FoldChange), probs = 0.999)
    fc.lower <- -fc.upper
    res <- res %>% 
        filter(baseMean <= mean.upper & log2FoldChange <= fc.upper & log2FoldChange >= fc.lower | gene %in% res.label$gene)

    # set filtering parameters
    if (ngene > 75000) {
        res.label <- res.label %>% 
            group_by(direction) %>% 
            slice_min(n = 50, order_by = padj, with_ties = F)}

    # set plot
    p <- res %>% 
        ggplot(aes(x = baseMean, y = log2FoldChange)) + 
        geom_point(aes(color = direction, size = size)) + 
        geom_text_repel(data = res.label, aes(label = gene), size = 3) + 
        geom_hline(yintercept = 0, linetype = "solid", size = 0.5) + 
        scale_size(range = c(0.5, 1.5)) + 
        scale_color_manual(values = c(Up = "red", Down = "blue", `Non-DE` = "grey60")) + 
        guides(color = guide_legend(title = "", override.aes = list(size = 5)), size = guide_none()) + 
        theme_border() + 
        theme_text() + 
        coord_cartesian(clip = "off") + 
        xlab("Log2 Mean Expression") + 
        ylab("Log2 Fold Change") + 
        ggtitle(paste0(comparison))

    # set caption
    caption <- paste0("total = ", as.character(label_comma()(ngene)), " features")
    p <- p + annotation_custom(grob = grid::textGrob(caption, 
        x = 1, y = -0.13, hjust = 1, gp = gpar(fontsize = 9, 
            col = "black")))

    # save plot
    if (save_plot) {
        save_plot(p, plot_name = paste0("diffexp_ma.pdf"), save_dir = save_dir, w = 7, h = 5)}
    print(p)
    return(p)
}

#' Generate Volcano Plot
#'
#' This function creates a volcano plot to visualize differential expression results.
#' It shows the relationship between statistical significance (-log10 adjusted p-value)
#' and biological significance (log2 fold change) for each gene.
#'
#' @param res Differential expression result data frame from run_diffexp()
#' @param n Number of top genes to label in each direction. Default is 25
#' @param fc.thresh Log2 fold change threshold for significance. Default is 1
#' @param p.thresh Adjusted p-value threshold for significance. Default is 0.05
#' @param crop Logical. If TRUE, limits the x-axis range to genes with significant changes. Default is TRUE
#' @param highlight.genes Vector of gene names to highlight. Default is NULL
#' @param save_plot Logical. If TRUE, saves the plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is current working directory
#' @return A ggplot object showing the volcano plot
#' @examples
#' \dontrun{
#' # Basic volcano plot
#' p <- plot_volcano(res)
#' 
#' # Volcano plot with custom thresholds and highlighted genes
#' p <- plot_volcano(res, n = 30, fc.thresh = 2, p.thresh = 0.01,
#'                   highlight.genes = c("GENE1", "GENE2"))
#' }
#' @export
plot_volcano <- function(res, n = 25, fc.thresh = 0.5, p.thresh = 0.05, crop = T, highlight.genes = NULL, save_plot = T, save_dir = getwd()){

    # validations
    res <- validate_res_comparison(res)
    validate_paths(save_dir)
    comparison <- unique(res$comparison)

    # set plot parameters
    line_cols <- "black"
    direction_cols <- c("red", "blue", "black")
    names(direction_cols) <- c("Up", "Down", "Non-DE")
    max.size <- 0.5
    ngene <- nrow(res)

    # add pxfc and rank columns
    res <- res %>%
        mutate(
            default = -log10(padj)*(log2FoldChange^2),
            direction = case_when(
                padj < p.thresh & log2FoldChange > fc.thresh ~ "Up",
                padj < p.thresh & log2FoldChange < -fc.thresh ~ "Down",
                .default = "Non-DE"),
            direction = factor(direction, c("Up", "Down", "Non-DE")))

    # crop x-axis to genes with significant changes
    if(crop){
        res.crop <- res %>% filter(padj < p.thresh)
        range <- sqrt(max(res.crop$log2FoldChange^2, na.rm = T))}
    else{
        range <- sqrt(max(res$log2FoldChange^2, na.rm = T))}

    # set plot parameters
    coord.y <- max(-log10(res$padj), na.rm = T)
    res <- res %>% 
        group_by(direction) %>%
        mutate(
            size = ifelse(direction == "Non-DE", 0.01, 0.5),
            count = n()*100/nrow(res),
            count = paste0(round(count, 1), "%"),
            count = ifelse(direction == "Non-DE", NA, count),
            coord.x = case_when(
                direction == "Up" ~ range*0.9,
                direction == "Down" ~ -range*0.9,
                .default = 0),
            coord.y = coord.y,
            ) %>%
        ungroup()

    # set highlight genes
    res.highlight <- res %>%
        filter(direction != "Non-DE") %>%
        group_by(direction) %>%
        slice_max(n = n, order_by = default)

    if(length(highlight.genes) > 0 & any(highlight.genes %in% res$gene)){
        max.size <- 3
        res <- res %>% 
            mutate(
                size = ifelse(gene %in% highlight.genes, 3, 0.01)) %>%
            ungroup()

        res.highlight <- res %>%
            filter(gene %in% highlight.genes)}

    # set plot
    p <- res %>%
        ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
        geom_vline(xintercept = 0, size = 0.4, color = "black") +
        geom_hline(yintercept = -log10(p.thresh), size = 0.4, color = line_cols, linetype = "dashed") +
        geom_vline(xintercept = -fc.thresh, size = 0.4, color = line_cols, linetype = "dashed") +
        geom_vline(xintercept = fc.thresh, size = 0.4, color = line_cols, linetype = "dashed") +
        geom_point(aes(color = direction, size = size)) +
        scale_size(range=c(0.01, max.size)) +
        geom_text(aes(x = coord.x, y = coord.y, label = count), size = 5) +
        scale_color_manual(values = direction_cols) +
        guides(color = guide_legend(title = "", override.aes = list(size = 5)), size = guide_none()) +
        geom_text_repel(data = res.highlight, aes(label = gene), size = 3.5) +
        xlim(c(-range, range)) +
        ggprism::theme_prism(border = T) +
        coord_cartesian(clip = 'off') +
        theme(plot.margin = margin(5,5,10,5)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        ggtitle(comparison)

    # set caption
    caption <- paste0("total = ", as.character(label_comma()(ngene)), " features")
    p <- p + annotation_custom(
        grob = grid::textGrob(caption, x = 1, y = -0.13, hjust = 1, gp = gpar(fontsize = 9, col = "black")))

    # save plot
    if(save_plot){
        save_plot(p, plot_name = paste0("diffexp_volcano.pdf"), save_dir = save_dir, w = 8, h = 5)}
    print(p)
    return(p)
}


#' Run Differential Expression Wrapper
#'
#' This function runs the differential expression analysis and generates MA and volcano plots for a single comparison
#'
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse".
#' @param group_by Column name in colData(dds) to use for grouping.
#' @param comparison Comparison to run.
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param one_to_all Logical. If TRUE, runs one-to-all comparisons. Default is FALSE
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return A data frame containing differential expression results for a comparison
#' @examples
#' \dontrun{
#' # Run differential expression pipeline
#' res <- run_diffexp_wrapper(dds)
#' 
#' # Run with custom grouping
#' res <- run_diffexp_wrapper(dds, group_by = "Treatment")
#' }
#' @export
run_diffexp_wrapper <- function(
    dds,
    org,
    group_by,
    comparison,
    order = "pxfc",
    save_dir = getwd()) {
    
    # validations
    dds <- validate_dds_comparison(dds, group_by, comparison)
    org <- validate_org(org)
    order <- validate_order(order)
    validate_paths(save_dir)

    # create save directory
    comparison_save_dir <- paste0(save_dir, "/", comparison)
    dir.create(comparison_save_dir, showWarnings = FALSE, recursive = TRUE)

    # Run differential expression analysis
    res <- run_diffexp(dds, org = org, group_by = group_by, comparison = comparison, order = order, save_dir = comparison_save_dir)

    # Generate MA plots
    plot <- plot_ma(res, save_plot = TRUE, save_dir = comparison_save_dir)
        
    # Generate volcano plots
    plot <- plot_volcano(res, save_plot = TRUE, save_dir = comparison_save_dir)

    return(res)}


#' Run Differential Expression Pipeline
#'
#' This function runs the differential expression analysis and generates MA and volcano plots for a single comparison
#'
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse".
#' @param group_by Column name in colData(dds) to use for grouping. 
#' @param one_to_all Logical. If TRUE, runs one-to-all comparisons. Default is FALSE
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param save_dir Directory where all output files will be saved. Default is current working directory
#' @return A data frame containing differential expression results for a comparison
#' @examples
#' \dontrun{
#' # Run differential expression pipeline
#' res <- run_diffexp_pip(dds, group_by = "group", one_to_all = T, save_dir = "results")
#' }
#' @export
run_diffexp_pip <- function(dds, org, group_by, order = "pxfc", one_to_all = FALSE, save_dir = getwd()){

    # validations   
    dds <- validate_dds_group_by(dds, group_by)
    org <- validate_org(org)
    order <- validate_order(order)
    validate_logical(one_to_all)
    validate_paths(save_dir)

    # initialize a list to store results
    res.list <- list()

    # run ONE-TO-ALL comparisons
    if(one_to_all & length(unique(dds[[group_by]])) > 2){

        message("Running ONE-TO-ALL differential expression analysis...")

        # create save directory
        group_save_dir <- paste0(save_dir, "/one-to-all_", group_by)
        dir.create(group_save_dir, showWarnings = FALSE, recursive = TRUE)

        # get group levels
        group.lv <- levels(dds[[group_by]])

        # run ONE-TO-ALL comparisons
        for(i in seq_along(group.lv)){
            dds.onetoall <- dds
            dds.onetoall[[group_by]] <- as.character(dds.onetoall[[group_by]])
            dds.onetoall[[group_by]] <- ifelse(dds.onetoall[[group_by]] == group.lv[i], group.lv[i], "All")
            dds.onetoall[[group_by]] <- factor(dds.onetoall[[group_by]], levels = c("All", group.lv[i]))
            design(dds.onetoall) <- design(dds)
            comparison_vec <- paste0(group.lv[i], "_vs_All")

            # run differential expression analysis
            res.list[[i]] <- run_diffexp_wrapper(dds.onetoall, org = org, group_by = group_by, comparison = comparison_vec, order = order, save_dir = group_save_dir)}}

    # run PAIRWISE comparisons
    else{
        message("Running PAIRWISE differential expression analysis...")

        # create save directory
        group_save_dir <- paste0(save_dir, "/pairwise_", group_by)
        dir.create(group_save_dir, showWarnings = FALSE, recursive = TRUE)

        # generate a vector of comparison
        group.lv <- levels(colData(dds)[[group_by]])
        comparison_vec <- generate_comparisons(group.lv)

        # run differential expression analysis pipeline for each comparison
        for(i in seq_along(comparison_vec)){

            # run differential expression analysis
            res.list[[i]] <- run_diffexp_wrapper(dds, org = org, group_by = group_by, comparison = comparison_vec[i], order = order, save_dir = group_save_dir)}}

    # merge results into a single data frame
    res <- bind_rows(res.list)
    return(res)
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