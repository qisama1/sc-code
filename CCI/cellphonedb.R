### cellphonedb

meta = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cellphonedb/mye-epi/meta.csv", row.names=1)

scRNA2 = scRNA[, rownames(meta)]

write.table(as.matrix(scRNA@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/9/cellphonedb/d.txt", sep='\t', quote=F)

write.table(as.matrix(cbind(rownames(scRNA@meta.data), scRNA@meta.data[,'ident', drop=F])), "/public/home/yuwenqi/sc-data/selected/9/cellphonedb/meta.txt", sep='\t', quote=F, row.names=F)

nohup python3 -u test.py >test.log 2>&1 &

nohup cellphonedb method statistical_analysis  meta.txt  d.txt     --counts-data=gene_name >test.log 2>&1 &
nohup cellphonedb method statistical_analysis  \
    meta.txt  \
    d.txt  \
    --counts-data gene_name  \
    --threshold 0.1 --iterations=10 --threads=8  >test.log 2>&1 & 


library('Matrix')
# Save normalised counts - NOT scaled!
writeMM(scRNA@assays$RNA@data, file = './endometrium_example_counts_mtx/matrix.mtx')
# save gene and cell names
write(x = rownames(scRNA@assays$RNA@data), file = "./endometrium_example_counts_mtx/features.tsv")
write(x = colnames(scRNA@assays$RNA@data), file = "./endometrium_example_counts_mtx/barcodes.tsv")

scRNA@meta.data$Cell = rownames(scRNA@meta.data)
df = scRNA@meta.data[, c('Cell', 'sub_clusters')]
write.table(df, file ='endometrium_example_meta.tsv', sep = '\t', quote = F, row.names = F)


nohup cellphonedb method statistical_analysis  \
    endometrium_example_meta.tsv  \
    endometrium_example_counts_mtx  \
    --counts-data gene_name  \
    --threshold 0.1 >test.log 2>&1 &
    
cellphonedb plot heatmap_plot "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/cellphonedb/endometrium_example_meta.tsv"