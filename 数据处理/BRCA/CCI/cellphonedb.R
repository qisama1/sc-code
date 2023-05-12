### cellphonedb
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", row.names=1)
scRNA = scRNA[, rownames(meta)]
scRNA[['ident']] = meta['sub_cluster']

module = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc_modules.csv")

scRNA_m1 = scRNA[, which(scRNA$ident %in% module$module1)]

write.table(as.matrix(scRNA_m1@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module1/d.txt", sep='\t', quote=F)

write.table(as.matrix(cbind(rownames(scRNA_m1@meta.data), scRNA_m1@meta.data[,'ident', drop=F])), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module1/meta.txt", sep='\t', quote=F, row.names=F)

scRNA_m1 = scRNA[, which(scRNA$ident %in% module$module2)]

write.table(as.matrix(scRNA_m1@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module2/d.txt", sep='\t', quote=F)

write.table(as.matrix(cbind(rownames(scRNA_m1@meta.data), scRNA_m1@meta.data[,'ident', drop=F])), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module2/meta.txt", sep='\t', quote=F, row.names=F)

scRNA_m1 = scRNA[, which(scRNA$ident %in% module$module3)]

write.table(as.matrix(scRNA_m1@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module3/d.txt", sep='\t', quote=F)

write.table(as.matrix(cbind(rownames(scRNA_m1@meta.data), scRNA_m1@meta.data[,'ident', drop=F])), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module3/meta.txt", sep='\t', quote=F, row.names=F)

scRNA_m1 = scRNA[, which(scRNA$ident %in% module$module4)]

write.table(as.matrix(scRNA_m1@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module4/d.txt", sep='\t', quote=F)

write.table(as.matrix(cbind(rownames(scRNA_m1@meta.data), scRNA_m1@meta.data[,'ident', drop=F])), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/module4/meta.txt", sep='\t', quote=F, row.names=F)


nohup cellphonedb method statistical_analysis  meta.txt  d.txt     --counts-data=gene_name >test.log 2>&1 &
nohup cellphonedb method statistical_analysis  \
    meta.txt  \
    d.txt  \
    --counts-data gene_name  \
    --threshold 0.1 --iterations=10 --threads=8  >test.log 2>&1 & 


library('Matrix')
setwd("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/cellphonedb/module1/matrix_method/")
# Save normalised counts - NOT scaled!
writeMM(scRNA_m1@assays$RNA@data, file = './mtx/matrix.mtx')
# save gene and cell names
write(x = rownames(scRNA@assays$RNA@data), file = "./mtx/features.tsv")
write(x = colnames(scRNA@assays$RNA@data), file = "./mtx/barcodes.tsv")

scRNA@meta.data$Cell = rownames(scRNA@meta.data)
df = scRNA@meta.data[, c('Cell', 'sub_clusters')]
write.table(df, file ='endometrium_example_meta.tsv', sep = '\t', quote = F, row.names = F)


nohup cellphonedb method statistical_analysis  \
    endometrium_example_meta.tsv  \
    endometrium_example_counts_mtx  \
    --counts-data gene_name  \
    --threshold 0.1 >test.log 2>&1 &
    
cellphonedb plot heatmap_plot "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/cellphonedb/endometrium_example_meta.tsv"