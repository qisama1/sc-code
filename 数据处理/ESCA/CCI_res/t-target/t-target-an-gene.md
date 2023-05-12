# needed_gene
```py
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ana_data.csv")
genes = set()
for gene in data.gene:
    genes.add(gene.split('|')[0])
    genes.add(gene.split('|')[1])
```

# get gene_df
```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/needed_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_df.csv")
```

# 处理python
```py
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_df.csv", index_col=0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['ident']
data['sample'] = meta['sample']

res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_gene_compare.csv", index_col = 0)
grade = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_gene_compare.csv", index_col=0)

gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ana_data.csv")
gene_sub = gene.loc[gene.type == 'ESCA']
for idx in gene_sub.index:
    v = gene_sub.loc[idx, ]
    gene_a = v.gene.split('|')[0]
    gene_b = v.gene.split('|')[1]
    cluster_1 = v.cluster.split('|')[0]
    cluster_2 = v.cluster.split('|')[1]
    p = pd.DataFrame(data.loc[data.sub_cluster == cluster_1].groupby('sample').mean().loc[:, gene_a])
    p.loc[:, 'type'] = res.loc[p.index, 'type']
    p_t = p.loc[p.type == 'Tumor']
    p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff/" + gene_a + "-" + cluster_1 + "_type.csv")
    p_t.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff/" + gene_a + "-" + cluster_1 + "_grade.csv")
    p = pd.DataFrame(data.loc[data.sub_cluster == cluster_2].groupby('sample').mean().loc[:, gene_b])
    p.loc[:, 'type'] = res.loc[p.index, 'type']
    p_t = p.loc[p.type == 'Tumor']
    p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff/" + gene_b + "-" + cluster_2 + "_type.csv")
    p_t.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff/" + gene_b + "-" + cluster_2 + "_grade.csv")

for idx in gene_sub.index:
    v = gene_sub.loc[idx, ]
    gene_a = v.gene.split('|')[0]
    gene_b = v.gene.split('|')[1]
    cluster_1 = v.cluster.split('|')[0]
    cluster_2 = v.cluster.split('|')[1]
    major = meta.loc[meta.cluster == meta.loc[data.loc[data.sub_cluster == cluster_1].index, 'cluster'].values[0]]
    p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, gene_a])
    p.loc[:, 'type'] = res.loc[p.index, 'type']
    p_t = p.loc[p.type == 'Tumor']
    p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff_major/" + gene_a + "-" + cluster_1 + "_type.csv")
    p_t.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff_major/" + gene_a + "-" + cluster_1 + "_grade.csv")
    major = meta.loc[meta.cluster == meta.loc[data.loc[data.sub_cluster == cluster_1].index, 'cluster'].values[0]]
    p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, gene_a])
    p.loc[:, 'type'] = res.loc[p.index, 'type']
    p_t = p.loc[p.type == 'Tumor']
    p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff_major/" + gene_b + "-" + cluster_2 + "_type.csv")
    p_t.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff_major/" + gene_b + "-" + cluster_2 + "_grade.csv")

```
# 画图
```R
setwd("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/gene_diff_major/")

fs = list.files('./', 'type')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('type'))
    boxplot_type(data, paste0("./plot/", filename, ".pdf"))
})

boxplot_type = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
    p=ggboxplot(data, x="variable", y="value", color = "type",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=type),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}

fs = list.files('./', 'grade')
lapply(fs, function(x) {
    print(x)
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('type', 'grade', 'stage'))
    boxplot_type2(data, paste0("./plot2/", filename, "_grade.pdf"))
})

lapply(fs, function(x) {
    print(x)
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('type', 'grade', 'stage'))
    boxplot_type3(data, paste0("./plot3/", filename, "_stage.pdf"))
})
boxplot_type2 = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"grade"])))]
    p=ggboxplot(data, x="variable", y="value", color = "grade",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=grade),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")), method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=grade),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}

boxplot_type3 = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"stage"])))]
    p=ggboxplot(data, x="variable", y="value", color = "stage",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=stage),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")), method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=stage),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}
```