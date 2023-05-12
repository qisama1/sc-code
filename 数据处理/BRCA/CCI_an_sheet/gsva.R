library(readr)
library(GSVA)
pathway_all <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_all2.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway_all <- as.data.frame(pathway_all)
pathway_list_all <- lapply(pathway_all, function(x) {
  unique(na.omit(x)) 
})

pathway_t <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_t.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway_t <- as.data.frame(pathway_t)
pathway_list_t <- lapply(pathway_t, function(x) {
  unique(na.omit(x)) 
})

pathway_mye <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_mye.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway_mye <- as.data.frame(pathway_mye)
pathway_list_mye <- lapply(pathway_mye, function(x) {
  unique(na.omit(x)) 
})

pathway_caf <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_caf.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway_caf <- as.data.frame(pathway_caf)
pathway_list_caf <- lapply(pathway_caf, function(x) {
  unique(na.omit(x)) 
})

pathway_end <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_end.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway_end <- as.data.frame(pathway_end)
pathway_list_end <- lapply(pathway_end, function(x) {
  unique(na.omit(x)) 
})

pathway_epi <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_epi.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway_epi <- as.data.frame(pathway_epi)
pathway_list_epi <- lapply(pathway_epi, function(x) {
  unique(na.omit(x)) 
})

scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", row.names=1)
scRNA = scRNA[, rownames(meta)]
scRNA[['ident']] =  meta$sub_cluster
Idents(scRNA) = scRNA$ident
epi = subset(scRNA, subset = cluster == 'Epi')
mye = subset(scRNA, subset = cluster == 'Mye')
tnk = subset(scRNA, subset = cluster == 'T/NK')
fib = subset(scRNA, subset = cluster == 'Fib')
end = subset(scRNA, subset = cluster == 'End')



ssgsea_score = gsva(as.matrix(as.data.frame(scRNA@assays$RNA@data)), pathway_list_all, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/all2.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(mye@assays$RNA@data)), pathway_list_mye, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/mye.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(tnk@assays$RNA@data)), pathway_list_t, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/t.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(epi@assays$RNA@data)), pathway_list_epi, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/epi.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(fib@assays$RNA@data)), pathway_list_caf, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/caf.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(end@assays$RNA@data)), pathway_list_end, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/end.csv")

