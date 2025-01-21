#时间轨迹分析，经典R包monocle3，基于基因表达，后续会继续退出PAGA时间轨迹分析软件
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
remove.packages('Matrix')
install.packages('irlba')
remotes::install_version("Matrix", version = "1.6-4")
install.packages("local/path/to/Matrix_1.6-4.tar.gz", repos = NULL, type="source")
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
#以上为R包安装部分，若已经安装可以忽略
setwd("G:/Ptbp1/crispr-mouse")

library(SeuratDisk)
#转换数据格式
Convert("G:/Ptbp1/crispr-mouse/adata_processed.nt.h5ad", "h5seurat",
        overwrite = TRUE,assay = "RNA")
scRNA <- LoadH5Seurat("G:/Ptbp1/crispr-mouse/adata_processed.nt.h5seurat")

meta <- scRNA@meta.data
# 使用 DimPlot 绘制 UMAP 图，按 cluster 分组


scRNA@meta.data$rename<-scRNA@meta.data$`Cluster-Name`
DimPlot(scRNA, reduction = "umap", group.by = "rename", label = TRUE, pt.size = 0.5) + NoLegend()
head(scRNA@meta.data)

pbmc <- scRNA
pbmc1 = pbmc [,pbmc @meta.data$rename %in% c("AT2-like","AT1-like","Gastric-like","High plasticity",
                                             "Lung progenitor-like","Endoderm-like","Early EMT-1","Pre-EMT")]
data <- GetAssayData(pbmc1, assay = 'RNA', slot = 'counts')

#其次获得meta信息
cell_metadata <- pbmc1@meta.data
head(cell_metadata)
#然后要获得gene_annotation 
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

#创建CDS对象
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#为啥要这样呢，这怎么解释呢，你用monocle的包，所以要听它的，还是要走一遍某些流程

#然后就是降维聚类分群
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 30)     #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
#像seurat一样展示pc数

#umap降维
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_pc_variance_explained(cds)   
#tSNE降维
cds <- reduce_dimension(cds, reduction_method="tSNE")
#cds <- reduce_dimension(cds, reduction_method="tSNE", perplexity = 10)

#聚类
#cds <- cluster_cells(cds) 
cds <- cluster_cells(cds, reduction_method = 'UMAP', cluster_method = "louvain")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc1, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed   #因此我们以后只要是画umap的点图就跟seurat的图片形状是一样的

#画图看一下
plot_cells(cds, reduction_method="UMAP", color_cells_by="rename")

cds <- learn_graph(cds, use_partition = T)


root_cells <- "L28.ATCTACTTCTAAGCCA-1" #指定root
cds <- order_cells(cds, root_cells = root_cells, reduction_method = "UMAP")
?order_cells
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds,
           color_cells_by = "rename",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=F,
           group_label_size=4,
           cell_size=1.5)
#以时间轨迹表示
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=FALSE,
           label_branch_points=F,
           cell_size=1.5) 
cds@metadata["RNA"]@pseudotime
cds@colData$pseudotime
pseudotime <- pseudotime[rownames(cds@colData)]
reducedDims <-cds@int_colData$reducedDims
cds@colData$pseudotime
summary(cds@colData$pseudotime)
pseudotime <- as.data.frame(pseudotime(cds,reduction_method = "UMAP"))  

#提取拟时分析结果返回seurat对象
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(pbmc1@meta.data)]
pbmc1@meta.data
pbmc1@meta.data$pseudotime <- pseudotime
p = FeaturePlot(pbmc1, reduction = "umap", features = "pseudotime")
# pseudotime中有无限值，无法绘图。
ggsave("Pseudotime_Seurat.pdf", plot = p, width = 8, height = 6)
saveRDS(pbmc1, file = "sco_pseudotime.rds")
write.csv(pbmc1@meta.data,file = "sco_pseudotime2.meta.csv",row.names = TRUE)
readRDS(pbmc1, file = "sco_pseudotime.rds")


#以下为不同格式的转换
saveRDS(pbmc1, file = "pbmc1.rds")
SaveH5Seurat(pbmc1, filename = "pbmc1.h5seurat", overwrite = TRUE)
Convert("pbmc1.h5seurat", "h5ad", overwrite = TRUE)
#monocle3至此结束