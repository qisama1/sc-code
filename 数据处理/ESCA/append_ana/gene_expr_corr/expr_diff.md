```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/an_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene_df.csv")
```

```py
# 1. 获取每行，然后判断类型，做一个简单的if/else即可，对于不同的类型对应不同的细胞类型的group
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['ident']
gene_used = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/an_gene.csv")
for i in gene_used.index:
    v = gene_used.iloc[i, ]
    g = v.gene
    if (g not in data.columns) :
        continue
    if (v.type == 'TAM'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('Macro')].index,]
        res = data.loc[major.index, [g, 'sub_cluster']]
        res.columns = ['gene', 'sub_cluster']
        res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene/" + v.type + "_" + g + ".csv")
    if (v.type == 'CD8+T'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('CD8')].index,]
        res = data.loc[major.index, [g, 'sub_cluster']]
        res.columns = ['gene', 'sub_cluster']
        res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene/" + v.type + "_" + g + ".csv")
    if (v.type == 'Epi'):
        major = meta.loc[meta.cluster == 'Epi']
        res = data.loc[major.index, [g, 'sub_cluster']]
        res.columns = ['gene', 'sub_cluster']
        res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene/" + v.type + "_" + g + ".csv")
    if (v.type == 'CAF'):
        major = meta.loc[meta.cluster == 'Fib']
        res = data.loc[major.index, [g, 'sub_cluster']]
        res.columns = ['gene', 'sub_cluster']
        res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene/" + v.type + "_" + g + ".csv")
    if (v.type == 'End'):
        major = meta.loc[meta.cluster == 'End']
        res = data.loc[major.index, [g, 'sub_cluster']]
        res.columns = ['gene', 'sub_cluster']
        res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene/" + v.type + "_" + g + ".csv")
```

```r
setwd("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene/")
fs = list.files('./', '.csv')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    marker = str_split(filename, '[_;]',simplify = T)[,2]
    data = read.csv(x, row.names=1)
    boxplot_diff(data, paste0("./plot/", filename, ".pdf"), marker)
})

boxplot_diff = function(data, filename, gene_name) {
    p1 = ggplot(data,aes(x=sub_cluster,y=gene),fill=sub_cluster)+
    scale_fill_manual(values=c('#ff890f','#7CFC00','#20B2AA',  '#F0F8FF',
                                '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                '#D9D9D9', '#EEAEEE'))+
    geom_boxplot(aes(fill=sub_cluster),show.legend = F)+ #outlier.shape = NA
    theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(),
            axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
            plot.title = element_text(hjust=0.5))+
    labs(x=NULL,
        y=NULL,
        title = gene_name)+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['gene']) - 0.15)  + stat_summary(fun='mean', geom="point")
    ggsave(filename, width = 1.8, height = 3, units = 'in')
}


```