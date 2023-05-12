
library(qs)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)

pbmc = qread("/public/home/yuwenqi/sc-data/selected/9/fib.qs")

setwd("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/monocle/")

##创建CDS对象并预处理数据
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
pbmc[['cell_type']] = pbmc[['ident']]

cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#cds$cell_type = meta[colnames(pbmc), 'sub_clusters']
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
cds <- cluster_cells(cds)
## 识别轨迹
cds <- learn_graph(cds, use_partition = FALSE)
pdf("./monocle_graph.pdf")
p = plot_cells(cds,
           color_cells_by = 'cell_type',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
p
dev.off()
## order cells
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,cds[['cell_type']] == 'MAFB+Epi']))

pdf('./order_cells.pdf')
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
dev.off()

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

pdf('./pseu_boxplot.pdf')
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(cell_type, monocle3_pseudotime, median), fill = cell_type)) +
  geom_boxplot()
dev.off()
## genes
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=4)


#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
                   pull(gene_short_name) %>% as.character()
#基因表达趋势图
write.csv(Track_genes, './track_genes-Basal.csv')
pdf("./genes_in_pseudotime-Basal.pdf")
plot_genes_in_pseudotime(cds[Track_genes_sig, ], color_cells_by="cell_type", 
                              min_expr=0.5, ncol = 2)
dev.off()