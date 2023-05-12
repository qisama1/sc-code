data = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/cellchat_key.csv", row.names=1)
scRNA = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_used.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all3.csv", row.names=1)
mye = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/mye.qs")
tnk = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/TNK/batch/cca/tnk.qs")
epi = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Epi/epi.qs")
fib = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Fib/fib.qs")
end = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/End/end.qs")

scRNA = scRNA[, rownames(meta)]
mye = mye[, rownames(meta)]
tnk = tnk[, rownames(meta)]
epi = epi[, rownames(meta)]
fib = fib[, rownames(meta)]
end = end[, rownames(meta)]

scRNA[['cluster']] = meta['cluster']
mye[['sub_cluster']] = meta['sub_cluster']
tnk[['sub_cluster']] = meta['sub_cluster']
fib[['sub_cluster']] = meta['sub_cluster']
epi[['sub_cluster']] = meta['sub_cluster']
end[['sub_cluster']] = meta['sub_cluster']

setwd("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_res/umap_cellchat/")

for (col in rownames(data)) {
    genes = strsplit(col, '[_;]')[[1]][1]
    gene1 = strsplit(genes, "[|;]")[[1]][1]
    gene2 = strsplit(genes, "[|;]")[[1]][2]
    f(gene1, gene2, str_replace(col,"[|;]","-"))
    if (length(strsplit(col, '[_;]')[[1]]) == 1) {
        next
    }
    gene3 = strsplit(col, '[_;]')[[1]][2]
    f(gene1, gene3, str_replace(col,"[|;]","-"))
}

f = function(gene1, gene2, filename) {
    path = paste0("./", filename)
    if (!dir.exists(path)) {
        dir.create(path)
    }
    p1 = FeaturePlot(scRNA, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(scRNA, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(scRNA, reduction = "umap", group.by = "cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_major.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(epi, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(epi, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(epi, reduction = "umap", group.by = "sub_cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Epi.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(end, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(end, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(end, reduction = "umap", group.by = "sub_cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_End.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(tnk, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(tnk, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(tnk, reduction = "umap", group.by = "sub_cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_TNK.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(mye, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(mye, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(mye, reduction = "umap", group.by = "sub_cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Mye.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(fib, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(fib, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(fib, reduction = "umap", group.by = "sub_cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Fib.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)


}
