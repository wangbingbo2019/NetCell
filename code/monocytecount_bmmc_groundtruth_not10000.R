library(caret)

#######################################
# 选取11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","09_pDC","10_cDC","17_B","19_CD8.N","20_CD4.N1","25_NK“细胞作为groundtruth
########################################

Seu_Healthy_Hema<-readRDS("E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")
meta.data<-Seu_Healthy_Hema@meta.data

celltypes<-c("11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","09_pDC","10_cDC",
             "17_B","19_CD8.N","20_CD4.N1","25_NK")

df<-meta.data[meta.data$BioClassification %in% celltypes,]
df$types<-"monocyte"
df$types[df$BioClassification=="09_pDC"]<-"non_monocyte"
df$types[df$BioClassification=="10_cDC"]<-"non_monocyte"
df$types[df$BioClassification=="19_CD8.N"]<-"non_monocyte"
df$types[df$BioClassification=="20_CD4.N1"]<-"non_monocyte"
df$types[df$BioClassification=="25_NK"]<-"non_monocyte"
df$types[df$BioClassification=="17_B"]<-"non_monocyte"

groundtruth_not1000<-Seu_Healthy_Hema[,rownames(df)]
groundtruth_not1000 <- FindVariableFeatures(groundtruth_not1000,nfeatures = 3000)
groundtruth_not1000 <- RunPCA(object = groundtruth_not1000, assay = "RNA", npcs = 50)

# 运行 UMAP
groundtruth_not1000_umap <- RunUMAP(groundtruth_not1000, dims = 1:50)
# 绘制 UMAP 图
png("E:/00_CeSOP/data/simulation/UMAP_plot_groundtruth_not1000.png", width = 800, height = 600)
DimPlot(groundtruth_not1000_umap, reduction = "umap", group.by = "BioClassification",pt.size=0.2,
        label = F, repel=TRUE)+scale_colour_manual(name = "BioClassification", values =color_scanpy_patient[c(2:9,1)])
dev.off()


# 可视化细胞得分
#lcczscore_file <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120.csv')
#lcczscore_file <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120-periphery.csv')
lcczscore_file <- read.csv('E:/00_CeSOP/results/simulation/Healthy-Hematopoiesis-191120-MAGMA.csv')

# 根据 list1 筛选 col1 中匹配的数据
filtered_df <- lcczscore_file[lcczscore_file$cell %in% rownames(df), ]
rownames(filtered_df)<-filtered_df$cell
# 按 groundtruth_umap 的行名重新排序 filtered_df
filtered_df <- filtered_df[rownames(groundtruth_not1000_umap@meta.data), ]
# 在 meta.data 中添加 lcczscore 列（示例值）
groundtruth_not1000_umap@meta.data$lcczscore <- filtered_df$lcczscore 
# 根据 lcczscore 连续值着色
dev.new()
png("E:/00_CeSOP/data/simulation/UMAP_plot_lcczscore_groundtruth_not1000.png", width = 800, height = 600)
FeaturePlot(groundtruth_not1000_umap, features = "lcczscore", pt.size = 0.1) +
  scale_color_gradientn(colors = c("white", "red"), 
                        limits = c(0, 10), name = "lcczscore")  
dev.off()



# 计算单核细胞占比
percent = sum(df$types == "monocyte")/nrow(df)
# 确定 lcczscore 的前 10% 阈值
threshold <- quantile(groundtruth_not1000_umap@meta.data$lcczscore, percent, na.rm = TRUE)
# 在 meta.data 中添加一个新列，标记是否为前 10%
groundtruth_not1000_umap@meta.data$top10 <- ifelse(groundtruth_not1000_umap@meta.data$lcczscore >= threshold, "Top10%", "Others")
# 提取 UMAP 坐标数据和元数据
umap_data <- Embeddings(groundtruth_not1000_umap, "umap") %>% as.data.frame()
umap_data$lcczscore <- groundtruth_not1000_umap@meta.data$lcczscore
umap_data$top10 <- groundtruth_not1000_umap@meta.data$top10


# 预测值——将 lcczscore 的前 10% 转化为二元值 1 和 0
groundtruth_not1000_umap@meta.data$pred <- ifelse(groundtruth_not1000_umap@meta.data$lcczscore >= threshold, 1, 0)
# 真实值——
groundtruth_not1000_umap@meta.data$true <- ifelse(df$types == "monocyte", 1, 0)

library(MLmetrics)
# 计算精确度
precision <- Precision(groundtruth_not1000_umap@meta.data$pred, groundtruth_not1000_umap@meta.data$true, positive = "1")
recall <- Recall(groundtruth_not1000_umap@meta.data$pred, groundtruth_not1000_umap@meta.data$true, positive = "1")
f1_score <- F1_Score(groundtruth_not1000_umap@meta.data$pred, groundtruth_not1000_umap@meta.data$true, positive = "1")
accuracy <- Accuracy(groundtruth_not1000_umap@meta.data$pred, groundtruth_not1000_umap@meta.data$true)

# 打印结果
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1 Score:", f1_score, "\n")
cat("Accuracy:", accuracy, "\n")


# 绘图：用透明度和颜色突出显示前 10% 的细胞
ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = lcczscore, alpha = top10), size = 1) +
  scale_color_gradientn(colors = c( "white", "red"), name = "LccZscore") +
  scale_alpha_manual(values = c("Top10%" = 1, "Others" = 0.5), guide = FALSE) +  # 设置透明度
  theme_minimal() +
  labs(title = paste("Precision:",round(precision, 3)))

