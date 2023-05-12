# 处理tnk
library(monocle)
library(qs)
library(Seurat)
setwd("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/monocle/")
pbmc = qread("/public/home/yuwenqi/sc-data/selected/9/tnk.qs")

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

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# 过滤低质量细胞
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

head(pData(HSMM))

# 选择基因
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~ident")
deg = subset(diff_test_res, qval < 0.01)
deg = deg[order(deg$qval, decreasing = F),]
ordergene = rownames(deg)[order(deg$qval)][1:200]

cds <- setOrderingFilter(HSMM, ordergene)

pdf("./train_ordergene.pdf")
plot_ordering_genes(cds)
dev.off()

cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds = orderCells(cds)