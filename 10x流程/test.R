
genes = read.csv()
data 

roc_gene = function(gene) {
    res<-roc(responses~BRCA[,gene],data=BRCA,aur=TRUE,
                     ci=TRUE, # 显示95%CI
                     #percent=TRUE, # 是否需要以百分比显示
                     smooth=TRUE,# 是否平滑曲线
                     levels=c('non_responders','responders'),direction=">" )
    return (res)
}

data = data.frame()

for (gene in b$gene_name) {
    res = roc_gene(gene)
    cur_data = data.frame(res$auc)
    rownames(cur_data) = gene
    data = rbind(data, cur_data)
}

