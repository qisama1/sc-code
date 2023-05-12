library(readr)
pathway <- read_delim("/public/home/yuwenqi/data/Data26/ppt23/genes.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway <- as.data.frame(pathway)
pathway_list <- lapply(pathway, function(x) {
  unique(na.omit(x)) 
})

exprSet_t = as.matrix(as.data.frame(scRNA@assays$RNA@data))
ssgsea_score = gsva(exprSet_t, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'

#write.csv(ssgsea_score, "/public/home/yuwenqi/data/Data26/ppt23/ssgsea/p19T_ssgsea.csv")
write.csv(ssgsea_score, "/public/home/yuwenqi/data/Data26/ppt23/ssgsea/p36T_ssgsea.csv")

library(GSVA)
library(GSEABase)
library(msigdbr)
library(limma)

## gsva
scRNA = readRDS("xx.rds")
exp = as.matrix(as.data.frame(scRNA@assays$RNA@data))

geneset <- getGmt('xxx.gmt')  
es <- gsva(as.matrix(exp), geneset, kcdf='Gaussian', method = 'gsva', parallel.sz=4, verbose = TRUE)

## ssgsea
exprSet_t = as.matrix(as.data.frame(scRNA@assays$RNA@data))
ssgsea_score = gsva(exprSet_t, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
