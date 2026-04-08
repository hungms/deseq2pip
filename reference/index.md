# Package index

## Pipeline Functions

- [`run_rna_pip()`](https://hungms.github.io/deseq2pip/reference/run_rna_pip.md)
  : Run Complete RNA-seq Pipeline
- [`run_atac_pip()`](https://hungms.github.io/deseq2pip/reference/run_atac_pip.md)
  : Run Complete ATAC-seq Pipeline
- [`deseq2pip-package`](https://hungms.github.io/deseq2pip/reference/deseq2pip-package.md)
  : DESeq2 Pipeline for RNA-seq Analysis
- [`enrichmentmap_pip()`](https://hungms.github.io/deseq2pip/reference/enrichmentmap_pip.md)
  : Format Enrichment Map Data
- [`plot_peak_annot_pip()`](https://hungms.github.io/deseq2pip/reference/plot_peak_annot_pip.md)
  : Plot ATAC-seq Peak Annotation
- [`run_diffexp_pip()`](https://hungms.github.io/deseq2pip/reference/run_diffexp_pip.md)
  : Run Differential Expression Wrapper
- [`run_dist_pip()`](https://hungms.github.io/deseq2pip/reference/run_dist_pip.md)
  : Run Sample Distance Pipeline This function runs the distance
  analysis portion of the RNA-seq analysis pipeline, including
  generating PCA and distance plots, and saving distance data.
- [`run_gsea_pip()`](https://hungms.github.io/deseq2pip/reference/run_gsea_pip.md)
  : Run GSEA pipeline
- [`run_qc_pip()`](https://hungms.github.io/deseq2pip/reference/run_qc_pip.md)
  : Run Quality Control Pipeline
- [`write_cls_pip()`](https://hungms.github.io/deseq2pip/reference/write_cls_pip.md)
  : Write class file pipeline for EnrichmentMap

## Wrapper Functions

## Quality Control (Module 1)

- [`remove_low_expression()`](https://hungms.github.io/deseq2pip/reference/remove_low_expression.md)
  : Remove Genes with Low Expression
- [`remove_mt_genes()`](https://hungms.github.io/deseq2pip/reference/remove_mt_genes.md)
  : Remove Mitochondrial Genes from DESeq2 Object
- [`remove_xy_genes()`](https://hungms.github.io/deseq2pip/reference/remove_xy_genes.md)
  : Remove XY Chromosome Genes from DESeq2 Object
- [`check_library()`](https://hungms.github.io/deseq2pip/reference/check_library.md)
  : Check Library Size Distribution
- [`run_batch_correction()`](https://hungms.github.io/deseq2pip/reference/run_batch_correction.md)
  : Run batch correction by Limma & Combat
- [`validate_batch()`](https://hungms.github.io/deseq2pip/reference/validate_batch.md)
  : Validate batch
- [`run_pca()`](https://hungms.github.io/deseq2pip/reference/run_pca.md)
  : Run Principal Component Analysis
- [`run_distance()`](https://hungms.github.io/deseq2pip/reference/run_distance.md)
  : Calculate and Plot Sample Distances

## Differential Expression (Module 2)

- [`read_diffexp()`](https://hungms.github.io/deseq2pip/reference/read_diffexp.md)
  : Read Differential Expression Results
- [`run_diffexp()`](https://hungms.github.io/deseq2pip/reference/run_diffexp.md)
  : Run Differential Expression Analysis
- [`run_diffexp_pip()`](https://hungms.github.io/deseq2pip/reference/run_diffexp_pip.md)
  : Run Differential Expression Wrapper

## Gene Set Enrichment Analysis (Module 3)

- [`plot_gsea_barplot()`](https://hungms.github.io/deseq2pip/reference/plot_gsea_barplot.md)
  : Plot GSEA Results
- [`read_gsea_rds()`](https://hungms.github.io/deseq2pip/reference/read_gsea_rds.md)
  : Read GSEA Results from RDS Files
- [`read_gsea_tsv()`](https://hungms.github.io/deseq2pip/reference/read_gsea_tsv.md)
  : Read GSEA Results from TSV Files
- [`run_gsea()`](https://hungms.github.io/deseq2pip/reference/run_gsea.md)
  : Run Gene Set Enrichment Analysis
- [`run_gsea_pip()`](https://hungms.github.io/deseq2pip/reference/run_gsea_pip.md)
  : Run GSEA pipeline
- [`validate_gsea_object()`](https://hungms.github.io/deseq2pip/reference/validate_gsea_object.md)
  : Validate gsea object
- [`validate_gsea_result()`](https://hungms.github.io/deseq2pip/reference/validate_gsea_result.md)
  : Validate gsea
- [`enrichmentmap_pip()`](https://hungms.github.io/deseq2pip/reference/enrichmentmap_pip.md)
  : Format Enrichment Map Data
- [`write_cls()`](https://hungms.github.io/deseq2pip/reference/write_cls.md)
  : Write class file pipeline for EnrichmentMap
- [`write_cls_pip()`](https://hungms.github.io/deseq2pip/reference/write_cls_pip.md)
  : Write class file pipeline for EnrichmentMap

## Plot Functions (Module 4)

- [`plot_gsea_barplot()`](https://hungms.github.io/deseq2pip/reference/plot_gsea_barplot.md)
  : Plot GSEA Results
- [`plot_ma()`](https://hungms.github.io/deseq2pip/reference/plot_ma.md)
  : Generate MA Plot
- [`plot_peak_annot()`](https://hungms.github.io/deseq2pip/reference/plot_peak_annot.md)
  : Plot ATAC-seq Peak Annotations
- [`plot_peak_annot_pip()`](https://hungms.github.io/deseq2pip/reference/plot_peak_annot_pip.md)
  : Plot ATAC-seq Peak Annotation
- [`plot_pie_chart()`](https://hungms.github.io/deseq2pip/reference/plot_pie_chart.md)
  : Create Pie Chart for ATACseq Annotation
- [`plot_volcano()`](https://hungms.github.io/deseq2pip/reference/plot_volcano.md)
  : Generate Volcano Plot

## Helper Functions

- [`import_msigdbr()`](https://hungms.github.io/deseq2pip/reference/import_msigdbr.md)
  : Import MSigDB Gene Sets
- [`import_nfcore_atac()`](https://hungms.github.io/deseq2pip/reference/import_nfcore_atac.md)
  : Import nfcore ATAC-seq DESeq2 Object
- [`import_nfcore_rna()`](https://hungms.github.io/deseq2pip/reference/import_nfcore_rna.md)
  : Import nfcore/rnaseq DESeq2 Object
- [`getTSS()`](https://hungms.github.io/deseq2pip/reference/getTSS.md) :
  Get TSS peaks from DESeq2 object
- [`summarize_genes_dds()`](https://hungms.github.io/deseq2pip/reference/summarize_genes_dds.md)
  : Summarize Gene Expression in DESeq2 Object
- [`generate_comparisons()`](https://hungms.github.io/deseq2pip/reference/generate_comparisons.md)
  : Generate Comparison Names
- [`validate_comparison()`](https://hungms.github.io/deseq2pip/reference/validate_comparison.md)
  : Validate DESeq2 object for comparison
- [`validate_comparisons()`](https://hungms.github.io/deseq2pip/reference/validate_comparisons.md)
  : Validate comparisons

## Validation Functions

- [`validate_atac()`](https://hungms.github.io/deseq2pip/reference/validate_atac.md)
  [`validate_dds_atac()`](https://hungms.github.io/deseq2pip/reference/validate_atac.md)
  : Validate DESeq2 object for ATAC-seq
- [`validate_batch()`](https://hungms.github.io/deseq2pip/reference/validate_batch.md)
  : Validate batch
- [`validate_comparison()`](https://hungms.github.io/deseq2pip/reference/validate_comparison.md)
  : Validate DESeq2 object for comparison
- [`validate_comparisons()`](https://hungms.github.io/deseq2pip/reference/validate_comparisons.md)
  : Validate comparisons
- [`validate_dds()`](https://hungms.github.io/deseq2pip/reference/validate_dds.md)
  : Validate DESeq2 object
- [`validate_design()`](https://hungms.github.io/deseq2pip/reference/validate_design.md)
  : Validate design formula
- [`validate_gene_set()`](https://hungms.github.io/deseq2pip/reference/validate_gene_set.md)
  : Validate gene set
- [`validate_gsea_object()`](https://hungms.github.io/deseq2pip/reference/validate_gsea_object.md)
  : Validate gsea object
- [`validate_gsea_result()`](https://hungms.github.io/deseq2pip/reference/validate_gsea_result.md)
  : Validate gsea
- [`validate_logical()`](https://hungms.github.io/deseq2pip/reference/validate_logical.md)
  : Validate logical
- [`validate_method()`](https://hungms.github.io/deseq2pip/reference/validate_method.md)
  : Validate methods
- [`validate_min_count()`](https://hungms.github.io/deseq2pip/reference/validate_min_count.md)
  : Validate min_count
- [`validate_msigdbr()`](https://hungms.github.io/deseq2pip/reference/validate_msigdbr.md)
  : Validate msigdbr
- [`validate_order()`](https://hungms.github.io/deseq2pip/reference/validate_order.md)
  : Validate order
- [`validate_org()`](https://hungms.github.io/deseq2pip/reference/validate_org.md)
  : Validate organism
- [`validate_pals()`](https://hungms.github.io/deseq2pip/reference/validate_pals.md)
  : Validate pals
- [`validate_paths()`](https://hungms.github.io/deseq2pip/reference/validate_paths.md)
  : Validate paths
- [`validate_res()`](https://hungms.github.io/deseq2pip/reference/validate_res.md)
  : Validate differential expression results
- [`validate_shape()`](https://hungms.github.io/deseq2pip/reference/validate_shape.md)
  : Validate shape
- [`validate_var()`](https://hungms.github.io/deseq2pip/reference/validate_var.md)
  : Validate DESeq2 object for RNA-seq

## Save Functions

- [`save_expression()`](https://hungms.github.io/deseq2pip/reference/save_expression.md)
  : Save Expression Data from DESeq2 Object
- [`save_plot()`](https://hungms.github.io/deseq2pip/reference/save_plot.md)
  : Save ggplot Object as PDF
- [`save_tsv()`](https://hungms.github.io/deseq2pip/reference/save_tsv.md)
  : Save Data Frame as TSV File

## Logging Functions

- [`log_output()`](https://hungms.github.io/deseq2pip/reference/log_output.md)
  : Log output
- [`log_renv()`](https://hungms.github.io/deseq2pip/reference/log_renv.md)
  : Log renv Snapshot
