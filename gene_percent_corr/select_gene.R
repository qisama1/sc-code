library(qs)
library(Seurat)

scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3_t.qs")

genes = c("COL1A1", "SDC1", "CXCL12", "CXCR4", "FN1", "THBS1", "THBS2", "GRN", "CXCL12", "TNFRSF1A", "SORT1", "CXCR4", "CXCL14", "COL1A1", 
"JAM2", "F11R", "HSPG2", "DAG1", "CXCL12", "CXCR4", "COL4A1", "JAG1", "CD46", "TNFRSF1A", "TNFRSF1B", "NECTIN2", "TIGIT", "ICAM2", "ITGAL", "ICAM1",
"COL4A1", "NECTIN2", "PECAM1", "CD38", "CD86", "CTLA4", "TNF", "VSIR", "CXCL9", "CXCR3", "TNFSF13B", "HLA-DPB1", "TFRC", "ICOS", "TNFRSF1B", "P2RY6",
"NAMPT", "CD86", "CD28", "MDK", "SDC1", "SDC2", "LRP1", "LAMB2", "ITGA2", "COL1A1", "ITGAV", "ITGA5", "CXCL12", "DPP4", "CXCR4", "CXCL14", "ITGB2", "ICAM1",
"TNFSF13B", "HLA-DPB1", "CD40", "EFNA1", "EPHA4", "EFNA3", "EFNA1", "EPHA2", "EFNA4", "NECTIN2", "SELP", "SELPLG", "PDGFB", "LRP1", "CCL14", "CCR1",
"CXCL12", "CXCR4", "CD44", "SELE")

genes = unique(genes)

expr = t(as.matrix(scRNA[genes,]@assays$RNA@data)) 