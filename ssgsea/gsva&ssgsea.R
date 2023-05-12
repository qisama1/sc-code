library(GSVA)
library(GSEABase)
library(msigdbr)
library(limma)
library(Seurat)
## gsva
scRNA = readRDS("xx.rds")
exp = as.matrix(as.data.frame(scRNA@assays$RNA@data))

geneset <- getGmt('xxx.gmt')  # gmt基因集读取
es <- gsva(as.matrix(exp), geneset, kcdf='Gaussian', method = 'gsva', parallel.sz=4, verbose = TRUE)

## ssgsea
exprSet_t = as.matrix(as.data.frame(scRNA@assays$RNA@data))
ssgsea_score = gsva(exprSet_t, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
