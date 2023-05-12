

data = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cellphonedb-key.csv", row.names=1)
for (col in rownames(data)) {
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

    p1 = FeaturePlot(Epi, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(Epi, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(Epi, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Epi.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(End, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(End, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(End, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_End.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(TNK, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(TNK, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(TNK, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_TNK.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(Mye, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(Mye, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(Mye, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Mye.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)

    p1 = FeaturePlot(Fib, reduction = "umap", features = c(gene1))
    p2 = FeaturePlot(Fib, reduction = "umap", features = c(gene2))
    p3 <- DimPlot(Fib, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)
    ggsave(paste0(path, "/", gene1 , "-", gene2, "_Fib.pdf"), plot=plot_grid(p1, p2, p3, nrow = 1), device='pdf', width=15, height=5)


}
