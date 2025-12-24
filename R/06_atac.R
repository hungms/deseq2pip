#' Plot ATAC-seq Peak Annotation
#'
#' This function plots the annotation distribution of ATAC-seq peaks.
#'
#' @param dds DESeq2 object containing the ATAC-seq peaks data
#' @param res Differential expression results data frame from run_diffexp()
#' @param save_dir Directory to save the plot. Default is current working directory
#' @param ... Additional arguments (ignored, for compatibility)
#' @return A list of ggplot objects showing the pie charts
#' @export
plot_peak_annot_pip <- function(
    dds,        
    res,
    save_dir = getwd(),
    ...) {
    
    # validations
    dds <- validate_dds_atac(dds)
    res <- validate_res(res)

    # get comparison
    comparisons <- unique(res$comparison)

    # run plot_peak_annot
    plots <- list()
    for(i in seq_along(comparisons)){
        selected_res <- res %>% filter(comparison == comparisons[i])
        plots[[i]] <- plot_peak_annot(dds, selected_res, save_dir = save_dir)}

    # return plots
    return(plots)
}

#' Get TSS peaks from DESeq2 object
#'
#' This function extracts TSS peaks from a DESeq2 object and removes duplicate genes.
#'
#' @param dds A DESeq2 object
#' @return A DESeq2 object containing only TSS peaks
#' @examples
#' \dontrun{
#' # Get TSS peaks
#' dds.tss <- getTSS(dds)
#' }
#' @export
getTSS <- function(dds){
    dds <- validate_dds_atac(dds)
    tss.dds <- dds[which(rowData(dds)$TSS == TRUE), ]
    duplicated.genes <- rowData(tss.dds)$gene[duplicated(rowData(tss.dds)$gene)]
    tss.dds <- tss.dds[-c(which(rowData(tss.dds)$gene %in% duplicated.genes)),]
    rownames(tss.dds) <- rowData(tss.dds)$gene
    tss.dds <- validate_dds_atac(tss.dds)
    return(tss.dds)}


#' Create Pie Chart for ATACseq Annotation
#' 
#' This function creates a pie chart to visualize the annotation distribution of ATACseq peaks.
#' 
#' @param annots A vector of peak annotations
#' @param comparison Comparison name
#' @param subtitle Subtitle for the plot
#' @return A ggplot object showing the pie chart
#' @export
plot_pie_chart <- function (annots, comparison, subtitle = "") {

    if(length(annots) == 0){
        return()}
    
    str_to_title_v2 <- function(input_string) {
        str_replace_all(input_string, "\\b[a-z]", function(x) str_to_title(x))
    }

    message("\t- generating pie charts for ", length(annots), " ", subtitle, "...")

    annot_counts <- table(annots[!is.na(annots)])

    annot_df <- data.frame(Annotation = str_to_title_v2(names(annot_counts)), 
        Count = as.numeric(annot_counts))

    annot_df <- annot_df %>% 
        dplyr::mutate(Percent = Count/sum(Count) * 100, Label = paste0(round(Percent, 1), "%")) %>% 
        arrange(Percent) %>% 
        dplyr::mutate(Annotation = factor(Annotation, levels = unique(Annotation)))

    p <- ggplot(annot_df, aes(x = "", y = Percent, fill = Annotation)) + 
        geom_bar(stat = "identity", width = 1) + 
        coord_polar("y", start = 0) + 
        no_gridlines() + 
        no_axis_text() + 
        xlab(NULL) + 
        ylab(NULL) + 
        geom_text(aes(label = Label), position = position_stack(vjust = 0.5)) + 
        ggtitle(comparison, subtitle = subtitle) + 
        theme(
            plot.background = element_blank(), 
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
            plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5))
    return(p)
}


#' Plot ATAC-seq Peak Annotations
#'
#' This function analyzes differentially expressed ATAC-seq peaks and creates pie charts
#' of their genomic annotations. It filters peaks by adjusted p-value and log2 fold change,
#' cleans up annotation labels, and generates three pie charts: one for all DE peaks,
#' one for upregulated peaks, and one for downregulated peaks.
#'
#' @param dds DESeq2 object containing the ATAC-seq peaks data
#' @param res Differential expression results data frame from run_diffexp()
#' @param p.thresh Adjusted p-value threshold for significance. Default is 0.05
#' @param fc.thresh Log2 fold change threshold for significance (absolute value). Default is 0.5
#' @param save_dir Directory to save the plot. Default is current working directory
#' @param save_plot Logical. If TRUE, saves the plots to PDF. Default is TRUE
#' @return A list of ggplot objects showing the pie charts
#' @examples
#' \dontrun{
#' # Generate annotation pie charts for ATAC-seq peaks
#' plots <- plot_peak_annot(dds, res)
#' 
#' # Use custom thresholds
#' plots <- plot_peak_annot(dds, res, p.thresh = 0.01, fc.thresh = 1)
#' }
#' @export
plot_peak_annot <- function(
    dds,
    res,
    p.thresh = 0.05,
    fc.thresh = 0.5,
    save_plot = TRUE,
    save_dir = getwd()) {
    
    # validations
    dds <- validate_dds_atac(dds)
    res <- validate_res(res)
    validate_logical(save_plot)

    # get comparison
    comparison <- unique(res$comparison)
    comparison_save_dir <- paste0(save_dir, "/", comparison, "/")

    # Filter differentially expressed peaks
    message("\t- extracting annotations for all differentially expressed peaks...")
    res <- res %>%
        filter(padj < p.thresh & abs(log2FoldChange) > fc.thresh)

    if (nrow(res) == 0) {
        message("\t- no peaks pass the filtering criteria for ", comparison)
        return()}

    peak_annots <- rowData(dds)$Annotation
    names(peak_annots) <- rownames(rowData(dds))
    
    # Remove text after " (" in all annotations
    peak_annots <- gsub(" \\(.*", "", peak_annots)

    # Get annotations for the DE peaks
    de_annots <- peak_annots[names(peak_annots) %in% res$peaks]

    # Split into up and down regulated
    up_idx <- res$log2FoldChange > 0
    up_annots <- peak_annots[names(peak_annots) %in% res$peaks[up_idx]]
    down_annots <- peak_annots[names(peak_annots) %in% res$peaks[!up_idx]]
    
    # Create the three pie charts
    all_plot <- plot_pie_chart(de_annots, subtitle = "DE Peaks", comparison = comparison)
    up_plot <- plot_pie_chart(up_annots, subtitle = "Upregulated DE Peaks", comparison = comparison)
    down_plot <- plot_pie_chart(down_annots, subtitle = "Downregulated DE Peaks", comparison = comparison)
    
    # Save plots if requested
    if(save_plot) {
        
        # Save each plot to its own PDF file using save_plot function
        save_plot(all_plot, 
                  plot_name = "peak_annotation_all.pdf", 
                  save_dir = comparison_save_dir, 
                  w = 6, h = 5)
        
        save_plot(up_plot, 
                  plot_name = "peak_annotation_up.pdf", 
                  save_dir = comparison_save_dir, 
                  w = 6, h = 5)
        
        save_plot(down_plot, 
                  plot_name = "peak_annotation_down.pdf", 
                  save_dir = comparison_save_dir, 
                  w = 6, h = 5)
        
    }
    print(all_plot)
    print(up_plot)
    print(down_plot)
    
    # Return the plots
    return(list(
        all = all_plot,
        up = up_plot,
        down = down_plot
    ))
}
