# Session Info

``` r
suppressPackageStartupMessages(library(deseq2pip))
```

    ## Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
    ## 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] deseq2pip_2.1.0             ggprism_1.0.7              
    ##  [3] renv_1.2.0                  enrichplot_1.30.5          
    ##  [5] ashr_2.2-63                 sva_3.58.0                 
    ##  [7] BiocParallel_1.44.0         genefilter_1.92.0          
    ##  [9] mgcv_1.9-4                  nlme_3.1-168               
    ## [11] labeling_0.4.3              limma_3.66.0               
    ## [13] forcats_1.0.1               scales_1.4.0               
    ## [15] ggpubr_0.6.3                fgsea_1.36.2               
    ## [17] msigdbr_26.1.0              clusterProfiler_4.18.4     
    ## [19] pheatmap_1.0.13             DESeq2_1.50.2              
    ## [21] SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ## [23] MatrixGenerics_1.22.0       matrixStats_1.5.0          
    ## [25] GenomicRanges_1.62.1        Seqinfo_1.0.0              
    ## [27] IRanges_2.44.0              S4Vectors_0.48.1           
    ## [29] BiocGenerics_0.56.0         generics_0.1.4             
    ## [31] plotr_0.1.3                 ggplotify_0.1.3            
    ## [33] RColorBrewer_1.1-3          wesanderson_0.3.7          
    ## [35] ggrepel_0.9.8               viridis_0.6.5              
    ## [37] viridisLite_0.4.3           cowplot_1.2.0              
    ## [39] patchwork_1.3.2             ggplot2_4.0.2              
    ## [41] strpip_0.1.4                biomaRt_2.66.2             
    ## [43] data.table_1.18.2.1         magrittr_2.0.5             
    ## [45] rlang_1.2.0                 stringr_1.6.0              
    ## [47] tibble_3.3.1                tidyr_1.3.2                
    ## [49] dplyr_1.2.1                
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] splines_4.5.3           filelock_1.0.3          R.oo_1.27.1            
    ##   [4] polyclip_1.10-7         XML_3.99-0.23           lifecycle_1.0.5        
    ##   [7] httr2_1.2.2             mixsqp_0.3-54           rstatix_0.7.3          
    ##  [10] edgeR_4.8.2             lattice_0.22-9          MASS_7.3-65            
    ##  [13] backports_1.5.1         sass_0.4.10             rmarkdown_2.31         
    ##  [16] jquerylib_0.1.4         yaml_2.3.12             ggtangle_0.1.1         
    ##  [19] DBI_1.3.0               abind_1.4-8             purrr_1.2.1            
    ##  [22] R.utils_2.13.0          yulab.utils_0.2.4       tweenr_2.0.3           
    ##  [25] rappdirs_0.3.4          gdtools_0.5.0           irlba_2.3.7            
    ##  [28] tidytree_0.4.7          annotate_1.88.0         pkgdown_2.2.0          
    ##  [31] codetools_0.2-20        DelayedArray_0.36.1     DOSE_4.4.0             
    ##  [34] ggforce_0.5.0           tidyselect_1.2.1        aplot_0.2.9            
    ##  [37] farver_2.1.2            BiocFileCache_3.0.0     jsonlite_2.0.0         
    ##  [40] Formula_1.2-5           survival_3.8-6          systemfonts_1.3.2      
    ##  [43] tools_4.5.3             ggnewscale_0.5.2        progress_1.2.3         
    ##  [46] treeio_1.34.0           ragg_1.5.2              Rcpp_1.1.1             
    ##  [49] glue_1.8.0              gridExtra_2.3           SparseArray_1.10.10    
    ##  [52] xfun_0.57               qvalue_2.42.0           withr_3.0.2            
    ##  [55] fastmap_1.2.0           truncnorm_1.0-9         digest_0.6.39          
    ##  [58] R6_2.6.1                gridGraphics_0.5-1      textshaping_1.0.5      
    ##  [61] GO.db_3.22.0            RSQLite_2.4.6           R.methodsS3_1.8.2      
    ##  [64] fontLiberation_0.1.0    prettyunits_1.2.0       httr_1.4.8             
    ##  [67] htmlwidgets_1.6.4       S4Arrays_1.10.1         scatterpie_0.2.6       
    ##  [70] pkgconfig_2.0.3         gtable_0.3.6            blob_1.3.0             
    ##  [73] S7_0.2.1                XVector_0.50.0          htmltools_0.5.9        
    ##  [76] fontBitstreamVera_0.1.1 carData_3.0-6           png_0.1-9              
    ##  [79] ggfun_0.2.0             knitr_1.51              reshape2_1.4.5         
    ##  [82] curl_7.0.0              cachem_1.1.0            parallel_4.5.3         
    ##  [85] AnnotationDbi_1.72.0    desc_1.4.3              pillar_1.11.1          
    ##  [88] vctrs_0.7.2             car_3.1-5               tidydr_0.0.6           
    ##  [91] dbplyr_2.5.2            xtable_1.8-8            cluster_2.1.8.2        
    ##  [94] evaluate_1.0.5          invgamma_1.2            cli_3.6.5              
    ##  [97] locfit_1.5-9.12         compiler_4.5.3          crayon_1.5.3           
    ## [100] SQUAREM_2026.1          ggsignif_0.6.4          plyr_1.8.9             
    ## [103] fs_2.0.1                ggiraph_0.9.6           stringi_1.8.7          
    ## [106] assertthat_0.2.1        babelgene_22.9          Biostrings_2.78.0      
    ## [109] lazyeval_0.2.3          GOSemSim_2.36.0         fontquiver_0.2.1       
    ## [112] Matrix_1.7-4            hms_1.1.4               bit64_4.6.0-1          
    ## [115] statmod_1.5.1           KEGGREST_1.50.0         igraph_2.2.3           
    ## [118] broom_1.0.12            memoise_2.0.1           bslib_0.10.0           
    ## [121] ggtree_4.0.5            fastmatch_1.1-8         bit_4.6.0              
    ## [124] ape_5.8-1               gson_0.1.0
