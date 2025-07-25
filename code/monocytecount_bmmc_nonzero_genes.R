library("scPagwas")
library("Seurat")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("stringr") 

# 读取单细胞表达数据，获取每个细胞的非零表达基因集合
scRNA_Healthy_Hema<-readRDS("E:/00_CeSOP/data/simulation/scRNA-Healthy-Hematopoiesis-191120.rds")
counts_matrix <- GetAssayData(scRNA_Healthy_Hema, slot = "counts")   # 获取 counts 矩阵
non_zero_genes_list <- apply(counts_matrix > 0, 2, function(x) rownames(counts_matrix)[x])  # 通过 apply() 函数获取每个细胞中非零表达的基因名称
non_zero_genes_df <- data.frame(
  Cell = names(non_zero_genes_list),     # 细胞名称
  Genes = sapply(non_zero_genes_list, paste, collapse = ", ") # 将每个基因列表合并为字符串
)
write.csv(non_zero_genes_df, "E:/00_CeSOP/data/simulation/Seu_Hema_non_zero_genes.csv", row.names = FALSE)


