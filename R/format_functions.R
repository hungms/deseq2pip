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