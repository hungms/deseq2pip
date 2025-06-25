pkgs <- c(
    "plotr", "DESeq2", "biomaRt", "ggalt", "pheatmap", "ggplotify", "clusterProfiler", "msigdbr", "msigdbdf", "fgsea", "ggpubr", "grid", "scales", "ggprism")

for(x in pkgs){
    usethis::use_package(x, type = "depends")}