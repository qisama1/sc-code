data1 = qread("/public/home/yuwenqi/sc-data/selected/9/tnk.qs")  #dataset10
data2 = qread("/public/home/yuwenqi/sc-data/selected/9/tnk.qs") #dataset9

meta1 = read.csv("/public/home/yuwenqi/sc-data/selected/9/tnk_info.csv", row.names=1)
meta2 = read.csv("/public/home/yuwenqi/sc-data/selected/9/tnk_info.csv", row.names=1)

data1$ident = meta1['ident']
data2$ident = meta2['ident']

ada1_df = data1$ident
ada2_df = data2$ident

gene1 <- VariableFeatures(data1)
gene2 <- VariableFeatures(data2)


ada1_exp <- as.matrix(data1@assays$RNA@data)
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- as.matrix(data2@assays$RNA@data)
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2))


res <- glm.predict(ada2_exp, ada2_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada1_exp, ada1_df, alpha = 0.99, nfolds = 10)

library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("#e9e9e9","white", "red"))
glm.predict.mean <-
      apply(res$logits, 2, function(e)
        sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
plot <- Heatmap(
      t(glm.predict.mean.prob),
	  col = col_fun,
      name = 'Predicted\nSimilarity',
      column_title = 'data10',
      row_title = 'data9',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
pdf("/public/home/yuwenqi/sc-data/selected/9-10-combined/9-9-combine_tnk.pdf")
plot
dev.off()