data = read.csv("xx")
min_cnt = 33
data = t(data)
sample = read.csv("xx")

gene_selected = data.frame()
for (gene in colnames(data)) {
    col_data = data[, gene]
    num = quantile(col_data, 0.4)
    col_data = col_data[sample$type == 'Tumor', ]
    cnt = sum(col_data >= num)
    if (cnt >= min_cnt) {
        t = data.frame(gene)
        rownames(t) = c(gene)
        colnames(t) = c('genes')
        gene_selected = rbind(gene_selected, t)
    }
}

