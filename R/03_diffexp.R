#' Run Differential Expression Wrapper
#'
#' This function runs the differential expression analysis and generates MA and volcano plots for a single comparison
#'
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse".
#' @param var Column name in colData(dds) to use for grouping.
#' @param comparisons Comparison to run.
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param one_to_all Logical. If TRUE, runs one-to-all comparisons. Default is FALSE
#' @param group_save_dir Directory where all output files will be saved. Default is current working directory
#' @return A data frame containing differential expression results for a comparison
#' @examples
#' \dontrun{
#' # Run differential expression pipeline
#' res <- run_diffexp_wrapper(dds)
#' 
#' # Run with custom grouping
#' res <- run_diffexp_wrapper(dds, var = "Treatment")
#' }
#' @export
run_diffexp_pip <- function(
    dds,
    design,
    org,
    var,
    comparisons = NULL,
    order = "pxfc",
    save_dir = getwd()) {
    
    # validations
    dds <- validate_var(dds, var)
    dds <- validate_design(dds, design)
    org <- validate_org(org)
    order <- validate_order(order)
    comparisons <- validate_comparisons(dds, var, comparisons)
    
    message("\n########################################################\nRunning Differential Expression Analysis Pipeline\n########################################################")
    message("running differential expression analysis for the following comparisons:\n")
    
    # Initialize list to store results
    res.list <- list()
    
    # run differential expression analysis for each comparison
    for(c in seq_along(comparisons)){
        
        # create save directory
        comparison_save_dir <- paste0(save_dir, "/", comparisons[c])
        dir.create(comparison_save_dir, showWarnings = FALSE, recursive = TRUE)

        # run differential expression analysis
        res <- run_diffexp(dds = dds, design = design, org = org, var = var, comparison = comparisons[c], order = order, save_data = TRUE, save_dir = comparison_save_dir)

        # generate ma and volcano plots
        plot_ma(res, save_plot = TRUE, save_dir = comparison_save_dir)
        plot_volcano(res, save_plot = TRUE, save_dir = comparison_save_dir)

        res.list[[c]] <- res}

    res.list <- bind_rows(res.list)
    return(res.list)}

#' Generate Comparison Names
#' 
#' This function generates comparison names for a given vector of variable levels.
#' It takes a vector of variable levels and returns a vector of comparison names.
#' 
#' @param var_levels A vector of variable levels
#' @return A vector of comparison names
#' @export
generate_comparisons <- function(var_levels){

    # generate pairwise comparisons
    var_levels <- unique(var_levels)

    if(length(var_levels) == 2){
        comparison_vec <- paste0(var_levels[2], "_vs_", var_levels[1])}
    else{
        combos <- combn(var_levels, 2)
        pairwise_vec <- apply(combos, 2, function(pair) {
            if (pair[2] < pair[1]) paste0(c(as.character(pair[2]), as.character(pair[1])), collapse = "_vs_")
            else paste0(c(as.character(pair)), collapse = "_vs_")})

        # generate one-to-all comparisons
        onetoall_vec <- c()
        for(i in seq_along(var_levels)){
            onetoall_vec[i] <- paste0(var_levels[i], "_vs_", paste0(var_levels[-i], collapse = "+"))}
        comparison_vec <- c(pairwise_vec, onetoall_vec)
            }
    
    # return comparison names
    return(comparison_vec)
}

#' Fail early with a clear message if the design matrix is rank-deficient.
#' Uses the same construction as DESeq2::designAndArgChecker when betaPrior is FALSE.
#' @noRd
assert_full_rank_model_matrix <- function(dds_obj, comparison, group_var) {
    des <- design(dds_obj)
    if (!inherits(des, "formula")) {
        return(invisible(NULL))
    }
    mm <- stats::model.matrix(des, data = as.data.frame(colData(dds_obj)))
    rk <- qr(mm)$rank
    nc <- ncol(mm)
    if (rk >= nc) {
        return(invisible(NULL))
    }
    pooled <- grepl("+", comparison, fixed = TRUE)
    msg <- paste0(
        "Model matrix is not full rank for comparison \"", comparison, "\" ",
        "(rank ", rk, " < ", nc, " columns)."
    )
    if (pooled) {
        msg <- paste0(
            msg,
            "\n\nOne-vs-rest comparisons merge several levels of `", group_var,
            "` into one group. Other design terms (e.g. batch) can then become ",
            "confounded with that merged grouping.\n\n",
            "Try: use a simpler design if appropriate (e.g. `~ ", group_var, "`); ",
            "or run only pairwise comparisons by passing `comparisons` without `+` ",
            "(e.g. `comps <- generate_comparisons(levels(colData(dds)[['",
            group_var, "']])); comps[!grepl(\"+\", comps, fixed = TRUE)]`)."
        )
    } else {
        msg <- paste0(
            msg,
            "\n\nCheck that terms in the design are not linearly redundant with `",
            group_var, "` for this contrast."
        )
    }
    stop(msg, call. = FALSE)
}

#' Run Differential Expression Analysis
#' 
#' This function runs the differential expression analysis for a single comparison.
#' 
#' @param dds DESeq2 object containing the expression data
#' @param org The organism to use, either "human" or "mouse".
#' @param var Column name in colData(dds) to use for grouping.
#' @param design Design formula for the DESeq2 object.
#' @param comparison Comparison to run.
#' @param order Column name to use for ranking genes. Default is "pxfc"
#' @param save_data Logical. If TRUE, saves the results to a TSV file. Default is TRUE
#' @param save_dir Directory to save the results. Default is current working directory
#' @return A data frame containing differential expression results for a comparison
#' @export
run_diffexp <- function(dds, org, var, design, comparison, order, save_data = TRUE, save_dir = getwd()){

    # validations
    dds <- validate_comparison(dds, var, comparison)
    dds <- validate_design(dds, design)
    org <- validate_org(org)
    order <- validate_order(order)
    validate_logical(save_data)

    org.to <- ifelse(org == "human", "mouse", "human")
    message("<<", comparison, ">>")

    # set contrast
    pair <- str_split(comparison, "_vs_") %>% unlist(.)
    contrast_vec <- c(var, pair[1], pair[2])

    # Independent colData copy: pooled one-vs-rest must not mutate `dds`, which
    # run_diffexp_pip() reuses for every comparison (would break later contrasts).
    temp_dds <- dds
    cd_df <- as.data.frame(colData(dds))
    if (!identical(rownames(cd_df), colnames(dds))) {
        rownames(cd_df) <- colnames(dds)
    }
    colData(temp_dds) <- S4Vectors::DataFrame(cd_df)

    # fix groupings (one-vs-rest: merge levels using exact matches, not regex gsub)
    if (str_detect(pair[1], "\\+") | str_detect(pair[2], "\\+")) {
        group1_vec <- str_split(pair[1], "\\+") %>% unlist(.) %>% str_trim(.)
        group2_vec <- str_split(pair[2], "\\+") %>% unlist(.) %>% str_trim(.)
        v <- as.character(colData(temp_dds)[[var]])
        if (length(intersect(group1_vec, group2_vec)) > 0) {
            stop(
                "Comparison \"", comparison, "\" has overlapping group names ",
                "between the two sides after splitting on '+'.",
                call. = FALSE
            )
        }
        v[v %in% group1_vec] <- pair[1]
        v[v %in% group2_vec] <- pair[2]
        if (!all(v %in% c(pair[1], pair[2]))) {
            leftover <- setdiff(unique(v), c(pair[1], pair[2]))
            stop(
                "After merging levels for comparison \"", comparison, "\", some samples ",
                "still have `", var, "` values not on either side of the contrast (e.g. ",
                paste(leftover, collapse = ", "), "). ",
                "Pooled comparisons must list every original level on one side or the other.",
                call. = FALSE
            )
        }
        colData(temp_dds)[[var]] <- factor(v, levels = c(pair[2], pair[1]))
    }

    # Match DESeq2: design slot must be set before model.matrix rank checks
    design(temp_dds) <- as.formula(design)
    assert_full_rank_model_matrix(temp_dds, comparison, var)

    # fit DESeq2 model
    message("\t- fitting DESeq2 model...")
    temp_dds <- quiet(DESeq(temp_dds))
    
    # extract differential expression results
    message("\t- extracting differential expression results...")
    res <- quiet(results(temp_dds, contrast = contrast_vec))

    # shrink lfc by ashr (contrast is already specified in res object)
    message("\t- shrinking lfc by ashr...")
    res <- quiet(lfcShrink(temp_dds, res = res, type = "ashr")) %>% 
        as.data.frame(.)    
    
    # add row metadata to results
    row.meta <- rowData(temp_dds) %>% as.data.frame(.)
    row.meta[which(colnames(row.meta) %in% colnames(res))] <- NULL
    if(ncol(row.meta) > 0){
        if(all(rownames(res) == rownames(assay(temp_dds)))){
            res <- res[rownames(assay(temp_dds)),]
            res <- cbind(res, row.meta)}

    # add pxfc and rank columns
    res <- res %>%
        mutate(
            comparison = comparison,
            pxfc = -log10(padj)*log2FoldChange,
            pxfc = ifelse(is.infinite(pxfc), 0, pxfc),
            rank = !!sym(order))

    # add gene annotations and filter results
    res <- quiet(run_annotation(res, org.from = org, org.to = org.to, gene_column = "gene")) %>%
        filter(!is.na(padj)) %>%
        filter(padj != "NA") %>%
        arrange(desc(rank))}

    # print number of DE genes
    all <- res %>% filter(padj < 0.05) %>% nrow()
    up <- res %>% filter(padj < 0.05 & log2FoldChange > 0) %>% nrow()
    down <- res %>% filter(padj < 0.05 & log2FoldChange < 0) %>% nrow()
    message(paste0("\t- ", all, " DE genes found: ", up, " upregulated and ", down, " downregulated"))

    if(save_data){
        save_tsv(res, tsv_name = paste0("diffexp_DESeq2.tsv"), save_dir = save_dir)}

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
    res <- validate_res(res)
    validate_logical(save_plot)
    comparison <- unique(res$comparison)
    message("\t- generating MA plot...")


    # set plot parameters
    quiet({
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
            var(direction) %>% 
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
    p <- quiet(p + annotation_custom(grob = grid::textGrob(caption, 
        x = 1, y = -0.13, hjust = 1, gp = gpar(fontsize = 9, 
            col = "black"))))

    # save plot
    if (save_plot) {
        save_plot(p, plot_name = paste0("diffexp_ma.pdf"), save_dir = save_dir, w = 7, h = 5)}
    print(p)
    })
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
    res <- validate_res(res)
    validate_logical(c(crop, save_plot))
    comparison <- unique(res$comparison)
    message("\t- generating volcano plot...")

    quiet({
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
    p <- quiet(p + annotation_custom(
        grob = grid::textGrob(caption, x = 1, y = -0.13, hjust = 1, gp = gpar(fontsize = 9, col = "black"))))

    # save plot
    if(save_plot){
        save_plot(p, plot_name = paste0("diffexp_volcano.pdf"), save_dir = save_dir, w = 8, h = 5)}
    print(p)
    })
    return(p)
}
