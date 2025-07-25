#if(!require(devtools)) install.packages("devtools")
#library(devtools)
#devtools::install_github("JSB-UCLA/scDesign2")
#install.packages("BiocManager")



library(scDesign2)
library(copula)    
library(Rtsne)
library(plyr)     
library(reshape2) 
library(gridExtra)
library(ggplot2); 
library(ggpubr)
library(cowplot)

library(stringr)
theme_set(theme_bw());

set.seed(42)
setwd("E:/00_CeSOP/data/simulation")

bulk_tpm<-read.delim("E:/00_CeSOP/data/simulation/pbmc_bulkdata/GSE107011_Processed_data_TPM.txt",header=T,row.names = 1)
Mono<-c("DZQV_C_mono","DZQV_I_mono","DZQV_NC_mono" ,"X925L_C_mono" ,"X925L_NC_mono","X9JD4_I_mono","X9JD4_NC_mono","G4YW_I_mono", "G4YW_NC_mono","X9JD4_mDC","G4YW_mDC")
dc<-c("DZQV_pDC","DZQV_mDC", "X925L_pDC","X925L_mDC" , "X9JD4_pDC","G4YW_pDC","G4YW_mDC","G4YW_I_mono", "G4YW_NC_mono")

B<- c("G4YW_B_NSM","G4YW_B_Ex","G4YW_B_SM","G4YW_B_naive","X9JD4_B_Ex","X9JD4_B_SM","X9JD4_B_NSM","X9JD4_B_naive","X925L_B_NSM","X925L_B_Ex","X925L_B_SM","X925L_B_naive","DZQV_B_SM","DZQV_B_naive","DZQV_B_NSM","DZQV_B_Ex")

NK<- c("DZQV_NK","X925L_NK","X9JD4_NK","G4YW_NK")

tcell<-c( "DZQV_CD8_naive","DZQV_CD8_CM","DZQV_CD8_EM","DZQV_CD8_TE","DZQV_TFH","DZQV_Treg","DZQV_Th1","DZQV_Th1.Th17","DZQV_Th17","DZQV_Th2","DZQV_CD4_naive", "X925L_CD8_naive","X925L_CD8_CM","X925L_CD8_EM","X925L_CD8_TE","X925L_TFH","X925L_Treg", "X925L_Th1","X925L_Th1.Th17","X925L_Th17", "X925L_Th2","X925L_CD4_naive","X925L_CD4_TE","X9JD4_CD8_naive", "X9JD4_CD8_CM","X9JD4_CD8_EM","X9JD4_CD8_TE" ,"X9JD4_TFH","X9JD4_Treg", "X9JD4_Th1","X9JD4_Th1.Th17","X9JD4_Th17","X9JD4_Th2","X9JD4_CD4_naive","X9JD4_CD4_TE","G4YW_CD8_naive","G4YW_CD8_CM" ,"G4YW_CD8_EM" ,"G4YW_CD8_TE","G4YW_TFH","G4YW_Treg" ,"G4YW_Th1", "G4YW_Th1.Th17","G4YW_Th17","G4YW_Th2","G4YW_CD4_naive","G4YW_B_naive","G4YW_B_NSM" )


bulk_tpm<-bulk_tpm[,c(Mono,dc,B,NK,tcell)]
genes<-unlist(str_sub(rownames(bulk_tpm),1,15))
bulk_tpm$gene<-genes
bulk_tpm<-bulk_tpm[!duplicated(bulk_tpm$gene),]
rownames(bulk_tpm)<-unique(bulk_tpm$gene)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")


library(org.Hs.eg.db)
library(dplyr)
a<-data.frame(ensembl_id=bulk_tpm$gene)
g2s=toTable(org.Hs.egSYMBOL)
g2e=toTable(org.Hs.egENSEMBL)
b=merge(a,g2e,by="ensembl_id",all.x=TRUE)
d=merge(b,g2s,by="gene_id",all.x=T)
d<-d[!duplicated(d$symbol),]
d<-d[!is.na(d$symbol),]

bulk_tpm<-bulk_tpm[d$ensembl_id,]
rownames(bulk_tpm)<-d$symbol
#save(bulk_tpm,file="/share/pub/dengcy/GWAS_Multiomics/pbmc_bulkdata/pbmc_bulk_tpm.RData")
bulk_tpm2<-bulk_tpm[,-ncol(bulk_tpm)]
bulk_tpm2<-ceiling(bulk_tpm2)
bulk_tpm2<-as.matrix(bulk_tpm2)
colnames(bulk_tpm2)<-c(rep("monocytes",length(Mono)),
                       rep("DC",length(dc)),
                       rep("B",length(B)),
                       rep("NK",length(NK)),
                       rep("Tcells",length(tcell)))
save(bulk_tpm2,file="E:/00_CeSOP/data/simulation/pbmc_bulkdata/pbmc_bulk_tpm3.RData")

library(SeuratObject)
library(Seurat)
library(scPagwas)
library(ggplot2)
load("E:/00_CeSOP/data/simulation/pbmc_bulkdata/pbmc_bulk_tpm3.RData")
bulk_tpm3<-CreateSeuratObject(counts=bulk_tpm2,assay = "RNA")

bulk_tpm3$celltype<-c(rep("monocytes",length(Mono)),
                      rep("DC",length(dc)),
                      rep("B",length(B)),
                      rep("NK",length(NK)),
                      rep("Tcells",length(tcell)))
bulk_tpm3 <- NormalizeData(bulk_tpm3, normalization.method = "LogNormalize", scale.factor = 10000)
bulk_tpm3 <- ScaleData(bulk_tpm3)
bulk_tpm3<-FindVariableFeatures(object=bulk_tpm3)
Idents(bulk_tpm3)<-bulk_tpm3$celltype



library(scDesign2)
library(copula)    
library(Rtsne)
library(plyr)     
library(reshape2) 
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggplot2); 
library(stringr)
theme_set(theme_bw());
library(scDesign2)

load("E:/00_CeSOP/data/simulation/pbmc_bulkdata/pbmc_bulk_tpm3.RData")

copula_result <- fit_model_scDesign2(bulk_tpm2,zp_cutoff=0.9,
                                     cell_type_sel=c('monocytes','DC','B','NK','Tcells'), 
                                     sim_method = 'copula',ncores = 1)
saveRDS(copula_result, file = 'E:/00_CeSOP/data/simulation/modelgroudtruth/copula_result_population2.rds')
sim_count_2000 <- simulate_count_scDesign2(copula_result, n_cell_new=2000, sim_method = 'copula',cell_type_prop = c(0.5,0.05,0.2,0.05,0.2))
rownames(sim_count_2000)<-rownames(bulk_tpm2)
saveRDS(sim_count_2000, file = 'sim_count_2000_population.rds')



#
library(SeuratObject)
library(Seurat)
sim_count_2000<-readRDS("E:/00_CeSOP/data/simulation/modelgroudtruth/sim_count_2000_population.rds")

sim_data<-CreateSeuratObject(counts=sim_count_2000,assay = "RNA")
sim_data$celltype<-colnames(sim_count_2000)
sim_data$type<- c(rep(1,1000),rep(0,1000))
sim_data <- NormalizeData(sim_data, normalization.method = "LogNormalize", scale.factor = 10000)
sim_data <- ScaleData(sim_data)
sim_data<-FindVariableFeatures(object=sim_data)

Idents(sim_data)<-sim_data$celltype
#sim_data<-RunTsne(object=sim_data)

sim_data <- FindVariableFeatures(sim_data,nfeatures = 3000)
sim_data <- RunPCA(object = sim_data, assay = "RNA", npcs = 50)

##subfunction:
cluster_pca_umap <- function(obj,assay=NULL, reduction,cluster_res = 0.3){
  #obj2 <- RunPCA(obj, assay = "SCT", reduction = "harmony",verbose = F)
  obj2 <- RunTSNE(object = obj,assay = assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- RunUMAP(object = obj2, assay =assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- FindNeighbors(object=obj2, assay = assay, reduction = reduction, dims = 1:50)
  obj2 <- FindClusters(object=obj2, resolution = cluster_res)
  return(obj2)
}
sim_data<-cluster_pca_umap(obj = sim_data,reduction="pca",cluster_res = 0.3)
saveRDS(sim_data,file="E:/00_CeSOP/data/simulation/modelgroudtruth/sim_data_8.16.rds")



