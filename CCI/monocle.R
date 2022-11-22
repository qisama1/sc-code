library(monocle)
library(Seurat)
library(dplyr)

pbmc = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3_t.qs")
expr_matrix = as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')

meta = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/bc3_meta_t_major2.csv", row.names=1)
pbmc[['cell_type']] = meta['sub_clusters']
p_data = pbmc@meta.data
p_data$celltype = pbmc@active.ident

f_data = data.frame(gene_short_name = row.names(pbmc), row.names = row.names(pbmc))

pd = new('AnnotatedDataFrame', data = p_data)
fd = new('AnnotatedDataFrame', data = f_data)

cds = newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())


# 估计size factor 离散度
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
qsave(cds, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/monocle/cds.qs")
express_genes = VariableFeatures(pbmc)
cds = setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

## 降维
cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')

## 模拟时间轴
cds = orderCells(cds)
pdf()
plot_cell_trajectory(cds, color_by = 'Pseudotime', size = 1, show_backbone = TRUE)
dev.off()

## 细胞类型上色
pdf()
plot_cell_trajectory(cds, color_by = "cell_type", size = 1, show_backbone = TRUE)


# expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
# cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
# gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))

# cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)


#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = pbmc@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#构建S4对象，CellDataSet
HSMM <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())