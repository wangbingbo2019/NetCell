library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer) 
set.seed(1)
epsilon = 1e-10
tissue_core = read.csv("E:\\00_CeSOP\\results\\asthma\\trait-cells-in-tissues-and-celltypes\\tissues\\top10_core_lcczscore.csv"
                       , na.strings = "", fileEncoding = "GB18030")
celltype_core = read.csv("E:\\00_CeSOP\\results\\asthma\\trait-cells-in-tissues-and-celltypes\\celltypes\\top10_core_lcczscore.csv"
                         , na.strings = "", fileEncoding = "GB18030")[,c('tissue','tissuename','celltype','celltypename','fisher.s.pvalue')]

# 找到 "column3" 小于 0.05 的行的索引
rows_to_change <- celltype_core$fisher.s.pvalue > 0.005
# 将符合条件的行的 "column2" 改为 "other"
celltype_core$celltype[rows_to_change] <- "其它"
celltype_core$fisher.s.pvalue[rows_to_change] <- 0.005
celltype_core = unique(celltype_core)

#d11 = data.frame(from="Asthma", to=tissue_core[['tissuename']])
d11 = data.frame(from="哮喘", to=tissue_core[['tissuename']])
d21 = data.frame(from=celltype_core[['tissuename']],
                 to=paste(celltype_core[['tissuename']],celltype_core[['celltypename']], sep = "---"))


edges<-rbind(d11[,1:2], d21[,1:2])

vertices_name<-unique(c(as.character(edges$from), as.character(edges$to)))
value1 <-rbind(1, tissue_core['fisher.s.pvalue'], celltype_core['fisher.s.pvalue'])
vertices<-data.frame(name = vertices_name, value = -log10(value1['fisher.s.pvalue'] + epsilon))
rownames(vertices)<-vertices_name

d2<-d21 %>% 
  mutate(order2=as.numeric(factor(from,
                                  levels=unique(from)[sort(summary (as.factor(from)),index.return=TRUE,decreasing = T)$ix],
                                  order=TRUE)))%>% 
  arrange(order2)

# 执行 left_join
d2 <- d2 %>%
  left_join(vertices, by = c("to" = "name"))
# 按照 order2 和 value 排序
d2 <- d2 %>%
  arrange(order2, desc(fisher.s.pvalue))

edges<-rbind(d11[,1:2], d2[,1:2])
list_unique<-unique(c(as.character(edges$from), as.character(edges$to)))
vertices = data.frame(
  name = list_unique, 
  value = vertices[list_unique,'fisher.s.pvalue']
) 

# 重新定义数据框，只保留 A 列值小于 2 的行
#vertices <- vertices %>% filter(value > -log10(0.01 + epsilon))
#edges <- edges[edges$from %in% vertices$name & edges$to %in% vertices$name, ]

vertices$group<-edges$from[match(vertices$name, edges$to)]

vertices$id<-NA
myleaves<-which(is.na( match(vertices$name, edges$from) ))
nleaves<-length(myleaves)
vertices$id<-nleaves/360*90
vertices$id[ myleaves]<-seq(1:nleaves)

vertices$angle<-90 - 360 * vertices$id / nleaves
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

mygraph <- graph_from_data_frame( edges, vertices, directed = TRUE )

# 设定内层颜色（例如红色）
#inner_layer_color <- "red"

# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_edge_diagonal(aes(colour=..index..)) +
#   scale_edge_colour_distiller(palette = "RdPu") +
#   geom_node_point(aes( x = x*1.07, y=y*1.07, fill=group, size=value*0.8), 
#                   shape=21,stroke=0.2,color='grey',alpha=0.8) +
#   geom_node_text(aes(x = x*1.25, y=y*1.25,  angle = angle,  label=name,color=group),
#                  size=2, alpha=1) +
#   scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
#   scale_fill_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
#   scale_size_continuous( range = c(0.1,7) ) +
#   expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))+
#   theme_void() +
#   theme(
#     legend.position="none",
#     plot.margin=unit(c(0,0,0,0),"cm"),
#     text = element_text(size = 11)  # 设置所有字体大小为11
#   )

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(aes(colour=..index..)) +
  scale_edge_colour_distiller(palette = "RdPu") +
  
  # 设定节点颜色（根据层级）
  geom_node_point(aes(x = x*1.07, y=y*1.07, 
                      fill = case_when(
                        name == "Asthma" ~ "black",       # 最内层设为黑色
                        name %in% d11$to ~ "grey",        # 第一层（组织层）设为灰色
                        TRUE ~ group                      # 第二层（细胞层）使用默认调色板
                      ), 
                      size=value*0.8), 
                  shape=21, stroke=0.2, color='grey', alpha=0.8) +
  
  # 设定文本颜色、字体大小和右对齐
  geom_node_text(aes(x = x*1.25, y=y*1.25, angle = angle, label=name, 
                     color = case_when(
                       name == "Asthma" ~ "black",   # 最内层文本黑色
                       name %in% d11$to ~ "grey",    # 第一层文本灰色
                       TRUE ~ group                  # 其余仍然使用默认调色板
                     )),
                 size=4, alpha=1) +  # 设置右对齐 (hjust = 1)
  
  # 调色板设置
  scale_fill_manual(values = c("black", "grey", rep(brewer.pal(20,"Paired"), 30))) +
  scale_colour_manual(values = c("black", "grey", rep(brewer.pal(20,"Paired"), 30))) +
  
  scale_size_continuous(range = c(0.1,7)) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
    text = element_text(size = 4)  # 设置所有字体大小为11
  )

ggsave("E:\\00_CeSOP\\plot\\figures\\asthma.png",
       width = 10,             # 宽
       height = 10,            # 高
       units = "in",          # 单位
       dpi = 500              # 分辨率DPI
)