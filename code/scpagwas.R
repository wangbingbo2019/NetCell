#######################################
# GWAS summary data progressed
########################################

library(readr)
setwd("E:/00_CeSOP/data/simulation")
############
monocytecount<-"E:/00_CeSOP/data/simulation/ieu-b-31.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table2(monocytecount,,comment = "#")

colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
  
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
})

value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
value_df<-as.data.frame(value_df)
colnames(value_df)<-c("ES","SE","LP","AF","SS")
gwas_data<-GWAS_raw[,c(1:5)]
gwas_data<- cbind(gwas_data,value_df)
gwas_data$LP <-  10^(-gwas_data$LP)
gwas_data$LP[which(gwas_data$LP>1)]<-1
colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
# if(is.na(gwas_data$maf[1])) 
# gwas_data$maf<-0.1
write.table(gwas_data,file="E:/00_CeSOP/data/simulation/monocytecount_gwas_data.txt",row.names=F,quote=F)


###################################################
# BMMC single cell data
###################################################
# library("scPagwas")
# library("Seurat")
# library("SingleCellExperiment")
# library("stringr") 
# scRNA_Healthy_Hema<-readRDS("E:/OneDrive/SingleCell/data/PBMCscATAC-seq/scRNA-Healthy-Hematopoiesis-191120.rds")
# counts <- assay(scRNA_Healthy_Hema, "counts")
# Seu_Healthy_Hema <- CreateSeuratObject(
#   counts = counts, 
#   meta.data=as.data.frame(colData(scRNA_Healthy_Hema)),
#   min.cells = 3, 
#   min.features = 200)
# 
# Idents(Seu_Healthy_Hema)<-scRNA_Healthy_Hema@colData$BioClassification
# table(Idents(Seu_Healthy_Hema))
# 
# #        01_HSC 02_Early.Eryth  03_Late.Eryth  04_Early.Baso    05_CMP.LMPP       06_CLP.1 
# #          1425           1653            446            111           2260            903 
# #        07_GMP    08_GMP.Neut         09_pDC         10_cDC 11_CD14.Mono.1 12_CD14.Mono.2 
# #          2097           1050            544            325           1800           4222 
# #  13_CD16.Mono         14_Unk       15_CLP.2       16_Pre.B           17_B      18_Plasma 
# #           292            520            377            710           1711             62 
# #      19_CD8.N      20_CD4.N1      21_CD4.N2       22_CD4.M      23_CD8.EM      24_CD8.CM 
# #          1521           2470           2364           3539            796           2080 
# #         25_NK         26_Unk 
# #          2143            161 
# 
# Seu_Healthy_Hema <- ScaleData(Seu_Healthy_Hema)
# Seu_Healthy_Hema <- NormalizeData(Seu_Healthy_Hema, normalization.method = "LogNormalize", scale.factor = 10000)
# saveRDS(Seu_Healthy_Hema,file="E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")



#############################################
# Run the single_data_input
#############################################
library(scPagwas)
Single_data<-readRDS("E:/00_CeSOP/data/simulation/Seu_Hema_data.rds")
Pagwas <- Single_data_input(Pagwas=NULL,
                            assay="RNA",
                            Single_data=Single_data,
                            Pathway_list=Genes_by_pathway_kegg)
Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                               Pathway_list=Genes_by_pathway_kegg
)
save(Pagwas,file="E:/00_CeSOP/data/simulation/Seu_Healthy_Hema_kegg_prePagwas.RData")



##############################################
# Run scPagwas
##############################################
library(scPagwas)
library(ggplot2)
suppressMessages(library(Seurat))
suppressMessages(library("dplyr"))
load("E:/00_CeSOP/data/simulation/Seu_Healthy_Hema_kegg_prePagwas.RData")
traits<-c("monocytecount")
# traits<-c("basophilcount","eosinophilcount" ,"Lymphocytecount3","monocytecount","neutrophilcount","WhiteBloodCellcount","LymphocytePercent","Hemoglobinconcen","MeanCorpuscularHemoglobin","MeanCorpusVolume")
for(i in traits){
  Pagwas<-scPagwas_main(Pagwas =Pagwas,
                        gwas_data =paste0("E:/00_CeSOP/data/simulation/",i,"_gwas_data.txt"),
                        Single_data ="E:/00_CeSOP/data/simulation/Seu_Hema_data_groundtruth5000.rds",
                        output.prefix=i,
                        Pathway_list=Genes_by_pathway_kegg,
                        output.dirs=paste0(i,"_scPagwas"),
                        assay="RNA",
                        block_annotation = block_annotation,
                        chrom_ld = chrom_ld)
  save(Pagwas,file=paste0("E:/00_CeSOP/data/simulation/",i,"_Hema_bmmc_scPagwas_groundtruth5000.RData"))
}


