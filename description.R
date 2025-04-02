pkgs <- c(
    "plotr", "DESeq2", "biomaRt", "ggalt", "pheatmap", "ggplotify", "clusterProfiler", "msigdbr", "msigdbdf", "fgsea", "ggpubr", "grid", "scales")

for(x in pkgs){
    usethis::use_package(x, type = "depends")}