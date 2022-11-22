pbmc = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/mye.qs")

setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/monocle/Mye/")

##创建CDS对象并预处理数据
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
pbmc[['cell_type']] = meta['sub_clusters']
pbmc[['cell_type']] = pbmc[['sub_cluster']]
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
pdf("./cds_umap.pdf")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") + ggtitle('cds.umap')
p1
dev.off()
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
pdf("./int_umap.pdf")
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") + ggtitle('int.umap')
p2
dev.off()
 
# monocle聚类 
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
        ggtitle("label by partitionID")
pdf("./monocle_cluster.pdf")
p = wrap_plots(p1, p2)
p
dev.off()

## 识别轨迹
cds <- learn_graph(cds)
pdf("./monocle_graph.pdf")
p = plot_cells(cds, color_cells_by = "cell_type",label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p
dev.off()

# 辅助线找root 
## cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(0,1,0.25))
embed <- data.frame(Embeddings(pbmc, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -3.8 & UMAP_1 < -3.5 & UMAP_2 >  -3.5 & UMAP_2 < -3)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
pdf("./pseudotime-Mono_1.pdf")
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
dev.off()