#' Plot Gene Expression
#'
#' This function creates boxplots showing the expression levels of selected genes across groups.
#' It can optionally add statistical significance indicators between groups.
#'
#' @param dds DESeq2 object containing the expression data
#' @param res Optional differential expression results data frame
#' @param group_by Column name in colData(dds) to group by.
#' @param genes Vector of gene names to plot
#' @param plot_name Name for the output plot
#' @param pal Vector of colors to use for groups. If NULL, uses default ggplot2 colors. Default is NULL
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
#' p <- plot_gene_exprs(dds, genes = c("GENE1", "GENE2"), plot_name = "my_genes", pal = c("red", "blue", "green"))
#' }
#' @export
plot_gene_exprs <- function(dds, res = NULL, group_by, genes, plot_name, pal = NULL, save_plot = T, save_dir = getwd()){
    
    #validations
    dds <- validate_dds_group_by(dds, group_by)
    res <- validate_res(res)
    pal <- validate_pal(pal)
    validate_paths(save_dir)

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

    if(!is.null(pal)) {
        p <- p + scale_fill_manual(values = pal)}

    if(!is.null(res)){
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

#' Plot GSEA Enriched Plot
#'
#' This function creates a plot to visualize GSEA results, showing the enriched gene set.
#' It can optionally add statistical significance indicators between groups.
#'
#' @param gsea GSEA result object from run_gsea()
#' @param gene_set Name of the gene set to plot
#' @param gene_set_title Title of the plot. Default is the gene set name
#' @param gene_set_title.size Size of the title text. Default is 8
#' @param show.pval Logical. If TRUE, shows p-value in the plot. Default is TRUE
#' @param show.fdr Logical. If TRUE, shows FDR in the plot. Default is TRUE
#' @param save_plot Logical. If TRUE, saves the plot to PDF. Default is TRUE
#' @param save_dir Directory to save the plot. Default is current working directory
#' @return A ggplot object showing the GSEA enrichment plot
#' @export
plot_gsea_enriched <- function(gsea.obj, gene_set, gene_set_title = NULL, gene_set_title_size = 10, show.pval = TRUE, show.fdr = TRUE, save_plot = T, save_dir = getwd()){

    #validations
    gsea.obj <- validate_gsea_object(gsea.obj)
    gene_set <- validate_gene_set(gsea.obj@result, gene_set)
    validate_paths(save_dir)

    id <- which(str_detect(gsea@result$ID, paste0("^", gene_set, "$")))
    plot <- gseaplot2(gsea, geneSetID = id, title = "", rel_heights = c(1, 0.2, 0.25))

    xmax <- min(plot[[1]]$data[[1]]) + (max(plot[[1]]$data[[1]]) - min(plot[[1]]$data[[1]]))*0.15
    ymax <- min(plot[[1]]$data[[2]]) + (max(plot[[1]]$data[[2]]) - min(plot[[1]]$data[[2]]))*0.2

    nes <- signif(gsea@result$NES[id], 3)
    label <- paste0("NES = ", nes)
    pval <- signif(gsea@result$pvalue[id], 1)
    fdr <- signif(gsea@result$qvalue[id], 1)

    if(length(gene_set_title) == 0){
        gene_set_title <- gene_set}
    
    if(nrow(gsea@result) < 30){
        show.fdr <- FALSE}
    if(show.pval){
        label <- paste0(label, "\np = ", pval)}
    if(show.fdr){
        label <- paste0(label, "\nFDR = ", fdr)}

    plot[[1]] <- plot[[1]] +
        geom_text(x = xmax, y = ymax, label = label, size = 3, fontface = "italic") +
        ylab("Enrichment Score") +
        theme_border() +
        theme_text() + 
        theme_gridlines() +
        ggtitle(
            comparison,
            subtitle = gene_set_title
            ) + 
        theme(
            panel.border = element_rect(size = 0.7, fill = NA, color = "black"),
            plot.subtitle = element_text(size = gene_set_title_size, hjust = 0.5, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position="none")
    
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
        theme(
            axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    
    plot <- as.ggplot(plot)
    if(save_plot){
        save_plot(plot, plot_name = paste0(comparison, "_", gene_set, ".pdf"), save_dir = paste0(save_dir, "/enrichplot/"), w = 5, h = 4)}

    return(plot)
}
