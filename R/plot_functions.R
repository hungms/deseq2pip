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
plot_ma <- function(res, save_plot = TRUE, save_dir = getwd()) {
  # Check if the input object meets the requirements
  stopifnot(is.data.frame(res) & all(c("baseMean", "log2FoldChange", "padj", "gene", "comparison") %in% colnames(res)))
  comparison <- unique(res$comparison)
  stopifnot(length(comparison) == 1)
  
  for(c in c("baseMean", "log2FoldChange", "padj")){
    res[[c]] <- as.numeric(res[[c]])}

  res <- res %>%
    mutate(baseMean = log2(baseMean + 1)) %>%
    filter(baseMean >= 3) %>%
    mutate(
        direction = case_when(
        log2FoldChange > 0 & padj < 0.05 ~ "Up",
        log2FoldChange < 0 & padj < 0.05 ~ "Down",
        .default = "Non-DE"),
      direction = factor(direction, c("Up", "Down", "Non-DE")),
      size = ifelse(direction == "Non-DE", 0.5, 1.5)
      )

  ngene <- nrow(res)
  res.label <- res %>%
    filter(direction != "Non-DE")

  # Remove outliers based on interquartile range (IQR)
  mean.upper <- quantile(res.label$baseMean, probs=0.999)
  fc.upper <- quantile(abs(res.label$log2FoldChange), probs=0.999)
  fc.lower <- -fc.upper
  
  res <- res %>%
    filter(
        baseMean <= mean.upper &
        log2FoldChange <= fc.upper &
        log2FoldChange >= fc.lower | gene %in% res.label$gene)
  
  if(ngene > 75000){
    res.label <- res.label %>%
        group_by(direction) %>%
        slice_min(n = 50, order_by = padj, with_ties = F)}

  caption <- paste0("total = ", as.character(label_comma()(ngene)), " features")

  # Create the ggplot
  p <- res %>%
    ggplot(aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = direction, size = size)) +
    geom_text_repel(data = res.label, aes(label = gene), size = 3) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5) +
    scale_size(range = c(0.5, 1.5)) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Non-DE" = "grey60")) +
    guides(color = guide_legend(title = "", override.aes = list(size = 5)), size = guide_none()) +
    theme_border() +
    theme_text() +
    coord_cartesian(clip = 'off') +
    xlab("Log2 Mean Expression") +
    ylab("Log2 Fold Change") +
    ggtitle(paste0("MA Plot; ", comparison))

  p <- p + annotation_custom(
    grob = grid::textGrob(caption, x = 1, y = -0.13, hjust = 1, gp = gpar(fontsize = 9, col = "black")))
  
  if(save_plot){
    save_plot(p, plot_name = paste0("diffexp_ma_plot.pdf"), save_dir = paste0(save_dir, "/", comparison), w = 7, h = 5)}
  
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
plot_volcano <- function(res, n = 25, fc.thresh = 1, p.thresh = 0.05, crop = T, highlight.genes = NULL, save_plot = T, save_dir = getwd()){
    comparison = unique(res$comparison)
    stopifnot(length(comparison) == 1)
    
    line_cols <- "black"
    direction_cols <- c("red", "blue", "black")
    names(direction_cols) <- c("Up", "Down", "Non-DE")
    max.size <- 0.5
    ngene <- nrow(res)

    if(max(abs(res %>% filter(padj < p.thresh) %>% .$log2FoldChange)) < fc.thresh){
        fc.thresh <- 0.5}

    res <- res %>%
        mutate(
            log2FoldChange = as.numeric(log2FoldChange),
            padj = as.numeric(padj),
            default = -log10(padj)*(log2FoldChange^2),
            direction = case_when(
                padj < p.thresh & log2FoldChange > fc.thresh ~ "Up",
                padj < p.thresh & log2FoldChange < -fc.thresh ~ "Down",
                .default = "Non-DE"),
            direction = factor(direction, c("Up", "Down", "Non-DE")))

    if(crop){
        res.crop <- res %>% filter(padj < p.thresh)
        range <- sqrt(max(res.crop$log2FoldChange^2, na.rm = T))}
    else{
        range <- sqrt(max(res$log2FoldChange^2, na.rm = T))}

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
            filter(gene %in% highlight.genes)
    }
    caption <- paste0("total = ", as.character(label_comma()(ngene)), " features")
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
        theme_border() +
        theme_text() +
        coord_cartesian(clip = 'off') +
        theme(plot.margin = margin(5,5,10,5)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        ggtitle(comparison)

    p <- p + annotation_custom(
        grob = grid::textGrob(caption, x = 1, y = -0.13, hjust = 1, gp = gpar(fontsize = 9, col = "black"))
        )

    print(p)

    if(save_plot){
        save_plot(p, plot_name = paste0("diffexp_volcano.pdf"), save_dir = paste0(save_dir, "/", comparison), w = 8, h = 5)
    }

    return(p)
}

#' Plot Gene Expression
#'
#' This function creates boxplots showing the expression levels of selected genes across groups.
#' It can optionally add statistical significance indicators between groups.
#'
#' @param dds DESeq2 object containing the expression data
#' @param res Optional differential expression results data frame
#' @param group_by Column name in colData(dds) to group by. Default is "Group1"
#' @param genes Vector of gene names to plot
#' @param plot_name Name for the output plot
#' @param cols Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
#' @param save_plot Logical. If TRUE, saves the plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is current working directory
#' @return A ggplot object showing gene expression boxplots
#' @examples
#' \dontrun{
#' # Basic expression plot
#' p <- plot_gene_exprs(dds, genes = c("GENE1", "GENE2"), plot_name = "my_genes")
#' 
#' # Expression plot with statistical significance
#' p <- plot_gene_exprs(dds, res = res, genes = c("GENE1", "GENE2"), plot_name = "my_genes")
#' 
#' # Expression plot with custom colors
#' p <- plot_gene_exprs(dds, genes = c("GENE1", "GENE2"), plot_name = "my_genes", cols = c("red", "blue", "green"))
#' }
#' @export
plot_gene_exprs <- function(dds, res = NULL, group_by = "Group1", genes, plot_name, cols = NULL, save_plot = T, save_dir = getwd()){

    stopifnot(is.factor(colData(dds)[[group_by]]))
    group.lv <- levels(colData(dds)[[group_by]])

    meta <- colData(dds) %>%
        as.data.frame(.) %>%
        mutate(samples = rownames(.)) %>%
        distinct(samples, !!sym(group_by))

    vsd <- vst(dds, blind = F)
    assay <- assay(vsd) %>%
        as.data.frame(.) %>%
        rownames_to_column("gene") %>%
        filter(gene %in% genes) %>%
        pivot_longer(!gene, names_to = "samples", values_to = "exprs") %>%
        merge(., meta, by = "samples", all.x = T)

    p <- assay %>%
        ggplot(aes_string(x = group_by, y = "exprs")) +
        geom_boxplot(aes_string(fill = group_by), width = 0.75) +
        geom_point() +
        facet_wrap(~gene, ncol = 5, scales = "free") +
        guides(fill = guide_legend(title = "")) +
        theme_border() +
        theme_text() +
        facet_aes() +
        xlab(NULL) +
        ylab("Expression") +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

    if(!is.null(cols)) {
        p <- p + scale_fill_manual(values = cols)
    }

    if(length(res) > 0){
        stopifnot(is.data.frame(res))
        yscale <- data.frame(group = group.lv, scale = seq(1, 1 + (length(group.lv)-1)*0.1, 0.1))

        ypos <- assay %>%
            filter(!!sym(group_by) %in% c(group.lv[1], group.lv[2])) %>%
            group_by(gene) %>%
            summarize(y.position = max(exprs)* 1.03) %>%
            ungroup()

        statdf <- res %>%
            filter(gene %in% genes) %>%
            mutate(
                `.y.` = "exprs",
                group1 = gsub("_vs_.*", "", comparison),
                group1 = factor(group1, group.lv),
                group2 = gsub(".*_vs_", "", comparison),
                group2 = factor(group2, group.lv),
                p = as.numeric(padj),
                label = case_when(
                    p < 0.001 ~ "***",
                    p < 0.01 ~ "**",
                    p < 0.05 ~ "*",
                    .default = "ns")
                ) %>%
            dplyr::select(c(".y.", "gene", "group1", "group2", "label")) %>%
            merge(., ypos, by = c("gene"), all.x = T) %>%
            group_by(gene) %>%
            arrange(gene, group2, group1) %>%
            mutate(
                rank = row_number(),
                y.position = y.position + (rank-1)*0.2) %>%
            ungroup() %>%
            dplyr::select(c(".y.", "gene", "group1", "group2", "label", "y.position"))

        ncomp <- choose(length(group.lv), 2)

        p <- p +
            stat_pvalue_manual(statdf, label = "label", tip.length = 0.015) +
            scale_y_continuous(expand = expansion(mult = c(0.1, ncomp*0.01)))
        }

    if(save_plot){
        ngene = length(genes)
        if(ngene > 5){
            ngenew = 5
            ngeneh = ceiling(ngene/5)}
        else{
            ngenew = ngene
            ngeneh = 1}
        ngroup = length(group.lv)

        w = (0.6*ngroup)*ngenew + 1
        h = 3.5*ngeneh

        if(length(res) > 0){
            h = (3.5 + 0.05*ncomp)*ngeneh}

        save_plot(p, plot_name = paste0(plot_name, "_exprs_vst.pdf"), save_dir = paste0(save_dir, "/plot/"), w = w, h = h)}

    print(p)
    return(p)
}

#' Plot GSEA Results
#'
#' This function creates a barplot to visualize GSEA results, showing the most significant
#' gene sets in each direction (up and down-regulated).
#'
#' @param gsea.df GSEA result data frame from run_gsea()
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
plot_gsea_barplot <- function(gsea.df, n = 10, signif = F, save_plot = T, save_dir = getwd()){
    comparison <- unique(gsea.df$comparison)
    collection <- unique(gsea.df$collection)
    stopifnot(length(comparison) == 1)

    stopifnot(all(c("NES", "qvalue", "ID") %in% colnames(gsea.df)))
    gsea.df$ID <- sub("_", "\\:", gsea.df$ID)
    gsea.df$NES <- as.numeric(gsea.df$NES)
    gsea.df$qvalue <- as.numeric(gsea.df$qvalue)

    selected_pathways <- gsea.df %>%
        arrange(desc(NES^2)) %>%
        mutate(
            direction = ifelse(NES > 0, "Up", "Down"),
            direction = factor(direction, c("Up", "Down"))) %>%
        group_by(direction) %>%
        slice_min(n = n, order_by = qvalue, with_ties = F) %>%
        .$ID

    if(signif){
        selected_pathways <- gsea.df %>%
            filter(qvalue < 0.05) %>%
            arrange(desc(NES^2)) %>%
            mutate(
                direction = ifelse(NES > 0, "Up", "Down"),
                direction = factor(direction, c("Up", "Down"))) %>%
            group_by(direction) %>%
            slice_min(n = n, order_by = qvalue, with_ties = F) %>%
            .$ID}

    selected.gsea.df <- gsea.df %>%
        filter(ID %in% selected_pathways) %>%
        select(ID, NES, pvalue, qvalue)

    empty_row <- data.frame(ID = "", NES = 0, pvalue = 1, qvalue = 1)

    yrange <- max(abs(selected.gsea.df$NES))*1.1
    
    p <- selected.gsea.df %>%
        rbind(., empty_row) %>%
        mutate(label = case_when(
            qvalue < 0.001 ~ "***",
            qvalue < 0.01 ~ "**",
            qvalue < 0.05 ~ "*",
            .default = "")) %>%
        mutate(
            size = ifelse(nchar(ID) >= 50, 1, ifelse(nchar(ID) >= 30, 2, 3)),
            ID =  str_replace_all(ID, ".{50}", "\\0\n")) %>%
        ggplot(aes(x = fct_reorder(ID, NES), y = NES, fill = NES)) +
        geom_col(aes(stroke = label), width = 0.75, col = "black") +
        geom_text(aes(y = ifelse(NES > 0, -0.1, 0.1), label = fct_reorder(ID, NES), hjust = ifelse(NES > 0, 1, 0), size = size), fontface = "bold") +
        scale_size(range = c(2, 3)) +
        geom_text(aes(label = label, y = NES + 0.2 * sign(NES)), position = position_dodge(width = 0.75), vjust = 0.75, size = 5) +
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
            size = guide_none()
            ) +
        ggtitle(paste0(collection, "; ", comparison)) +
        theme_border() +
        theme_text() +
        theme(
            panel.border = element_rect(fill = NA, color = "grey40", size = 0.3),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
        theme_gridlines() +
        ylim(c(-yrange, yrange)) +
        geom_hline(yintercept = 0, color = "black", size = 0.7) +
        scale_x_discrete(expand=c(0.05, 0.05)) +
        coord_flip()

    if(save_plot){
        plot_name <- paste0("gsea_", collection, "_top", n)
        if(signif){
            plot_name <- paste0(plot_name, "_signif")}
        save_plot(p, plot_name = paste0(plot_name, "_barplot.pdf"), save_dir = paste0(save_dir, "/", comparison, "/"), w = 9, h = 7)}
    print(p)
    return(p)
}

#' Plot GSEA Enriched Plot
#'
#' This function creates a plot to visualize GSEA results, showing the enriched gene set.
#' It can optionally add statistical significance indicators between groups.
#'
#' @param gsea GSEA result object from run_gsea()
#' @param gene_set Name of the gene set to plot
#' @param title.size Size of the title text. Default is 8
#' @param show.pval Logical. If TRUE, shows p-value in the plot. Default is TRUE
#' @param show.fdr Logical. If TRUE, shows FDR in the plot. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is current working directory
#' @return A ggplot object showing the GSEA enrichment plot
#' @export
plot_gsea_enriched <- function(gsea, gene_set, title.size = 8, show.pval = TRUE, show.fdr = TRUE, save_plot = T, save_dir = getwd()){
    
    comparison <- unique(gsea@result$comparison)
    stopifnot(length(comparison) == 1)
    stopifnot(gene_set %in% gsea@result$ID)

    id <- which(str_detect(gsea@result$ID, gene_set))
    plot <- gseaplot2(gsea, geneSetID = id, title = "", rel_heights = c(1, 0.2, 0.25))

    xmax <- min(plot[[1]]$data[[1]]) + (max(plot[[1]]$data[[1]]) - min(plot[[1]]$data[[1]]))*0.12
    ymax1 <- min(plot[[1]]$data[[2]]) + (max(plot[[1]]$data[[2]]) - min(plot[[1]]$data[[2]]))*0.2
    ymax2 <- min(plot[[1]]$data[[2]]) + (max(plot[[1]]$data[[2]]) - min(plot[[1]]$data[[2]]))*0.08

    nes <- signif(gsea@result$NES[id], 3)
    label <- paste0("NES = ", nes)
    pval <- signif(gsea@result$pvalue[id], 1)
    fdr <- signif(gsea@result$qvalue[id], 1)

    if(length(title) == 0){
        title <- gsea@result$ID[id]}
    
    if(nrow(gsea@result) < 100){
        show.fdr <- FALSE}
    if(show.pval){
        label <- paste0(label, "\np = ", pval)}
    if(show.fdr){
        label <- paste0(label, "\nFDR = ", fdr)}

    plot[[1]] <- plot[[1]] +
        geom_text(x = xmax, y = ymax1, label = comparison, size = 4, fontface = "bold") +
        geom_text(x = xmax, y = ymax2, label = label, size = 3) +
		geom_hline(yintercept = 0, linetype = "dashed") +
        ylab("Enrichment Score") +
        theme_border() +
        theme_text() + 
        theme_gridlines() +
        ggtitle(gene_set) + 
        theme(legend.position="none")
    
    plot[[2]] <- plot[[2]] +
        theme_border() +
        no_gridlines() +
        theme_text() +
        no_axis_text() + 
        theme(legend.position="none")
    
    plot[[3]] <- plot[[3]] +
        ylab("Rank") +
        theme_border() +
        no_gridlines() +
        theme_text() +
        no_axis_text()
    
    plot <- as.ggplot(plot)
    print(plot)
    if(save_plot){
        save_plot(plot, plot_name = paste0(comparison, "_", gene_set, ".pdf"), save_dir = paste0(save_dir, "/enrichplot/"), w = 5, h = 4)}

    return(plot)
}