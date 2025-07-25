library("scPagwas")
library("Seurat")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("stringr") 
library(ggplot2)

#读取单细胞表达数据，进行标准化处理
scRNA_Healthy_Hema<-readRDS("E:/00_CeSOP/data/simulation/scRNA-Healthy-Hematopoiesis-191120.rds")
counts <- assay(scRNA_Healthy_Hema, "counts")
Seu_Healthy_Hema <- CreateSeuratObject(
  counts = counts, 
  meta.data=as.data.frame(colData(scRNA_Healthy_Hema)),
  min.cells = 3, 
  min.features = 200)
Idents(Seu_Healthy_Hema)<-scRNA_Healthy_Hema@colData$BioClassification
table(Idents(Seu_Healthy_Hema))
Seu_Healthy_Hema <- NormalizeData(Seu_Healthy_Hema, normalization.method = "LogNormalize", scale.factor = 10000)
Seu_Healthy_Hema <- ScaleData(Seu_Healthy_Hema)
saveRDS(Seu_Healthy_Hema,file="E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")


#提取细胞名称和 BioClassification 
Seu_Healthy_Hema<-readRDS("E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")
cell_names <- rownames(Seu_Healthy_Hema@meta.data)
bio_classification <- Seu_Healthy_Hema@meta.data$BioClassification
cell_classification_df <- data.frame(Cell = cell_names, BioClassification = bio_classification)
write.csv(cell_classification_df, file = "E:/00_CeSOP/data/simulation/Seu_Hema_cell_classification.csv", row.names = FALSE, quote = FALSE)


# 标准化和 PCA 降维
Seu_Healthy_Hema<-readRDS("E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")
Seu_Healthy_Hema <- FindVariableFeatures(Seu_Healthy_Hema,nfeatures = 3000)
Seu_Healthy_Hema <- RunPCA(object = Seu_Healthy_Hema, assay = "RNA", npcs = 50)

# 绘制 t-SNE 图
# tsne <- RunTSNE(Seu_Healthy_Hema, dims = 1:50)
# png("E:/00_CeSOP/data/simulation/TSNE_plot_allcells.png", width = 800, height = 600)
# DimPlot(tsne, reduction = "tsne", group.by = "BioClassification")
# dev.off()

# 绘制 UMAP 图
umap <- RunUMAP(Seu_Healthy_Hema, dims = 1:50)
png("E:/00_CeSOP/data/simulation/UMAP_plot_allcells.png", width = 800, height = 600)
DimPlot(umap, reduction = "umap", group.by = "BioClassification",pt.size=0.2)
dev.off()


# 可视化细胞得分
#lcczscore_df <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120.csv')
#lcczscore_df <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120-periphery.csv')
#lcczscore_df <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120-MAGMA.csv')
lcczscore_df <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120-SCPAGWAS.csv')

rownames(lcczscore_df)<-lcczscore_df$cell
# 按 umap 的行名重新排序 filtered_df
lcczscore_df <- lcczscore_df[rownames(umap@meta.data), ]
# 在 meta.data 中添加 lcczscore 列（示例值）
umap@meta.data$lcczscore <- lcczscore_df$lcczscore 
# 根据 lcczscore 连续值着色
png("E:/00_CeSOP/data/simulation/UMAP_plot_lcczscore_allcells.png", width = 800, height = 600)
 
FeaturePlot(umap, features = "lcczscore", pt.size = 0.1) +
  scale_color_gradientn(colors = c("white", "red"), 
                        name = "lcczscore")  
dev.off()


#### 
meta.data<-Seu_Healthy_Hema@meta.data
df<-meta.data
df$types<-"non_monocyte"
df$types[df$BioClassification=="11_CD14.Mono.1"]<-"monocyte"
df$types[df$BioClassification=="12_CD14.Mono.2"]<-"monocyte"
df$types[df$BioClassification=="13_CD16.Mono"]<-"monocyte"


# 计算单核细胞占比
percent = sum(df$types == "monocyte")/nrow(df)
# 确定 lcczscore 的前 10% 阈值
threshold <- quantile(umap@meta.data$lcczscore, percent, na.rm = TRUE)
# 在 meta.data 中添加一个新列，标记是否为前 10%
umap@meta.data$top10 <- ifelse(umap@meta.data$lcczscore >= threshold, "Top10%", "Others")
# 提取 UMAP 坐标数据和元数据
umap_data <- Embeddings(umap, "umap") %>% as.data.frame()
umap_data$lcczscore <- umap@meta.data$lcczscore
umap_data$top10 <- umap@meta.data$top10


# 预测值——将 lcczscore 的前 10% 转化为二元值 1 和 0
umap@meta.data$pred <- ifelse(umap@meta.data$lcczscore >= threshold, 1, 0)
# 真实值——
umap@meta.data$true <- ifelse(df$types == "monocyte", 1, 0)

library(MLmetrics)
# 计算精确度
precision <- Precision(umap@meta.data$pred, umap@meta.data$true, positive = "1")
recall <- Recall(umap@meta.data$pred, umap@meta.data$true, positive = "1")
f1_score <- F1_Score(umap@meta.data$pred, umap@meta.data$true, positive = "1")
accuracy <- Accuracy(umap@meta.data$pred, umap@meta.data$true)
# 打印结果
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1 Score:", f1_score, "\n")
cat("Accuracy:", accuracy, "\n")

# 示例数据
y_true <- umap@meta.data$true  # 实际值
y_pred <- umap@meta.data$pred  # 预测值

# 去除 NA 值
valid_idx <- complete.cases(y_true, y_pred)
y_true <- y_true[valid_idx]
y_pred <- y_pred[valid_idx]

# 转换为因子
y_true <- factor(y_true, levels = c(0, 1))
y_pred <- factor(y_pred, levels = c(0, 1))

# 使用 confusionMatrix
library(caret)
conf_matrix <- confusionMatrix(data = y_pred, reference = y_true, positive = "1")

# 打印结果
print(conf_matrix)




# 
# # 最大最小归一化
# umap_data$lcczscore_normalized <- (umap_data$lcczscore - min(umap_data$lcczscore)) / 
#   (max(umap_data$lcczscore) - min(umap_data$lcczscore))
# # 按照 [0, 0.2, 0.4, 0.6, 0.8, 1] 分段
# umap_data$lcczscore_category <- cut(
#   umap_data$lcczscore_normalized,
#   breaks = c(0, 0.1, 0.3, 0.5, 1),
#   labels = c("0-0.1", "0.1-0.2", "0.2-0.4", "0.4--1"),
#   include.lowest = TRUE
# )
# # 定义颜色
# category_colors <- c("#fcfaf1", "#f6e7ce", "#f4c7a8", "#d67573")
# # 映射分段到颜色
# umap_data$lcczscore_color <- category_colors[as.numeric(umap_data$lcczscore_category)]
# # 绘制 UMAP 图，分段着色
# ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
#   geom_point(aes(color = lcczscore_color), size = 0.5) +
#   scale_color_identity() +  # 使用定义好的颜色
#   theme_minimal() +
#   labs(title = paste("Precision:",precision))


# 绘图：用透明度和颜色突出显示前 10% 的细胞
ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = lcczscore, alpha = top10), size = 0.2, alpha = 0.3) +
  scale_color_gradientn(colors = c( "white", "red"), name = "LccZscore") +
  scale_alpha_manual(values = c("Top10%" = 1, "Others" = 0), guide = FALSE) +  # 设置透明度
  theme_minimal() +
  labs(title = paste("Precision:",precision))
