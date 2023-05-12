# 读基因
```R
gene = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/needed_gene.csv")
write.csv(as.matrix(scRNA[gene$Gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/gene_df.csv")
```

```python
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/gene_df.csv").T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['sub_cluster']
data['sample'] = meta['orig.ident']

cc = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cellchat-key.csv", index_col=0)
res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_gene_scores/s_celltype_gene_type.csv", index_col = 0)
grade = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_gene_scores/s_celltype_gene_grade.csv", index_col=0)

for idx in cc.index:
    for col in cc.columns:
      if (cc.loc[idx, col] > 0):
        gene_a = idx.split('|')[0]
        gene_b = idx.split('|')[1]
        if (gene_b == 'ITGA5_ITGB1'):
            continue
        cluster_1 = col.split('|')[0]
        cluster_2 = col.split('|')[1]
        p = pd.DataFrame(data.loc[data.sub_cluster == cluster_1].groupby('sample').mean().loc[:, gene_a])
        p.loc[:, 'type'] = res.loc[p.index, 'type']
        p_t = p.loc[p.type == 'tumor']
        p_t.loc[:, 'grade'] = grade.loc[p_t.index, 'grade']
        p.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/boxplot_cellchat/" + gene_a + "-" + cluster_1 + "_type.csv")
        p_t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/boxplot_cellchat/" + gene_a + "-" + cluster_1 + "_grade.csv")
        p = pd.DataFrame(data.loc[data.sub_cluster == cluster_2].groupby('sample').mean().loc[:, gene_b])
        p.loc[:, 'type'] = res.loc[p.index, 'type']
        p_t = p.loc[p.type == 'tumor']
        p_t.loc[:, 'grade'] = grade.loc[p_t.index, 'grade']
        p.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/boxplot_cellchat/" + gene_b + "-" + cluster_2 + "_type.csv")
        p_t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/boxplot_cellchat/" + gene_b + "-" + cluster_2 + "_grade.csv")

```

# 读取所有的文件
```R
setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/boxplot_cellchat/")

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
    data=melt(data, id.vars=c('type', 'grade'))
    boxplot_type2(data, paste0("./plot2/", filename, "_grade.pdf"))
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
```