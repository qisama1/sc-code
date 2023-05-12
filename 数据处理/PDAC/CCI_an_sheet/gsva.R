library(readr)
library(GSVA)
pathway_all <- read_delim("/public/home/yuwenqi/sc-data/selected/pathway/pathway_all.csv", ",", 
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

scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", row.names=1)
mye = qread("/public/home/yuwenqi/sc-data/selected/35/mye.qs")
tnk = qread("/public/home/yuwenqi/sc-data/selected/35/tnk.qs")
epi = qread("/public/home/yuwenqi/sc-data/selected/35/epi.qs")
fib = qread("/public/home/yuwenqi/sc-data/selected/35/fib.qs")
end = qread("/public/home/yuwenqi/sc-data/selected/35/end.qs")
b = qread("/public/home/yuwenqi/sc-data/selected/35/bcell.qs")

scRNA = scRNA[, rownames(meta)]
mye = mye[, rownames(meta)]
tnk = tnk[, rownames(meta)]
epi = epi[, rownames(meta)]
fib = fib[, rownames(meta)]
end = end[, rownames(meta)]
b = b[, rownames(meta)]


ssgsea_score = gsva(as.matrix(as.data.frame(scRNA@assays$RNA@data)), pathway_list_all, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/all.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(mye@assays$RNA@data)), pathway_list_mye, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/mye.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(tnk@assays$RNA@data)), pathway_list_t, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/t.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(epi@assays$RNA@data)), pathway_list_epi, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/epi.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(fib@assays$RNA@data)), pathway_list_caf, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/caf.csv")

ssgsea_score = gsva(as.matrix(as.data.frame(end@assays$RNA@data)), pathway_list_end, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(ssgsea_score, "/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/end.csv")

