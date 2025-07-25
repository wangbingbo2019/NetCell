library("scPagwas")
library("Seurat")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("stringr") 
library(ggplot2)

############################
# 1. lcczscore,核心基因，
############################

# 1.1 单细胞数据可视化
Seu_Healthy_Hema<-readRDS("E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")
Seu_Healthy_Hema <- FindVariableFeatures(Seu_Healthy_Hema,nfeatures = 3000)
Seu_Healthy_Hema <- RunPCA(object = Seu_Healthy_Hema, assay = "RNA", npcs = 50)
umap <- RunUMAP(Seu_Healthy_Hema, dims = 1:50)
#png("E:/00_CeSOP/data/simulation/UMAP_plot_allcells.png", width = 800, height = 600)
DimPlot(umap, reduction = "umap", group.by = "BioClassification",pt.size=0.2)
#dev.off()


# 1.2 lcczscore可视化
lcczscore_df <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120.csv')
rownames(lcczscore_df)<-lcczscore_df$cell
lcczscore_df <- lcczscore_df[rownames(umap@meta.data), ]
umap@meta.data$lcczscore <- lcczscore_df$lcczscore 
# 根据 lcczscore 连续值着色
png("E:/00_CeSOP/data/simulation/UMAP_plot_lcczscore_allcells.png", width = 800, height = 600)

FeaturePlot(umap, features = "lcczscore", pt.size = 0.1) +
  scale_color_gradientn(colors = c("white", "red"), 
                        name = "lcczscore")  
dev.off()


df<-Seu_Healthy_Hema@meta.data
df$types<-"non_monocyte"
df$types[df$BioClassification=="11_CD14.Mono.1"]<-"monocyte"
df$types[df$BioClassification=="12_CD14.Mono.2"]<-"monocyte"
df$types[df$BioClassification=="13_CD16.Mono"]<-"monocyte"


# 计算单核细胞占比
percent = sum(df$types == "monocyte")/nrow(df)
# 确定 lcczscore 的前 10% 阈值
threshold <- quantile(umap@meta.data$lcczscore, 1-percent, na.rm = TRUE)

# 预测值——将 lcczscore 的前 10% 转化为二元值 1 和 0
pred <- ifelse(umap@meta.data$lcczscore >= threshold, 1, 0)
# 真实值——
true <- ifelse(df$types == "monocyte", 1, 0)


library(MLmetrics)
# 评估模型性能
conf_matrix <- confusionMatrix(factor(pred), factor(true))
# 获取并打印各项指标
accuracy <- conf_matrix$overall["Accuracy"]          # 准确率
precision <- conf_matrix$byClass["Pos Pred Value"]   # 精确率
recall <- conf_matrix$byClass["Sensitivity"]         # 召回率
f1_score <- conf_matrix$byClass["F1"]                # F1分数
specificity <- conf_matrix$byClass["Specificity"]    # 特异性
sensitivity <- conf_matrix$byClass["Sensitivity"]    # 灵敏度（与召回率相同）
fpr <- 1 - specificity                              # 假阳性率

# 打印各项指标
cat("Precision:", precision, "\n")
cat("Accuracy:", accuracy, "\n")
cat("False Positive Rate:", fpr, "\n")
