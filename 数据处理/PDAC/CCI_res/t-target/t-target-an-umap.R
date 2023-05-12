
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

scRNA[['cluster']] = meta['cluster']
mye[['sub_cluster']] = meta['ident']
tnk[['sub_cluster']] = meta['ident']
fib[['sub_cluster']] = meta['ident']
epi[['sub_cluster']] = meta['ident']
end[['sub_cluster']] = meta['ident']
b[['sub_cluster']] = meta['ident']

setwd("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/PDAC/umap/")

data = read.csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ana_data.csv")
data = subset(data, type == 'PDAC')
for (col in data$gene) {
    gene1 = strsplit(col, "[|;]")[[1]][1]
    gene2 = strsplit(col, "[|;]")[[1]][2]
    f(gene1, gene2)
}

f = function(gene1, gene2) {
    path = paste0("./", gene1, "-", gene2)
    if (!dir.exists(path)) {
        dir.create(path)
    }
    p1 = FeaturePlot(scRNA, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(scRNA, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(scRNA, reduction = "umap", group.by = "cluster",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_major.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(epi, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(epi, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(epi, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Epi.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(end, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(end, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(end, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_End.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(tnk, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(tnk, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(tnk, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_TNK.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(mye, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(mye, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(mye, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Mye.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(fib, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(fib, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(fib, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Fib.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)


}
