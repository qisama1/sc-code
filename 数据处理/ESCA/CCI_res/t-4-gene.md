# 读基因
```R
gene = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/needed_gene.csv")
write.csv(as.matrix(scRNA[gene$Gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/9/module/gene_df.csv")

```

```python
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/gene_df.csv", index_col=0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['ident']
data['sample'] = meta['sample']

cellchat_clu1 = ['Ecm_myCAF_1', 'Wound_myCAF', 'Angiogenic_1', 'Angiogenic_2', 'End_progenitor', 'Macro_SPP1']
cellchat_clu2 = ['NDRG1+Epi', 'Cycling_Epi', 'CD8Teff_2', 'CD8Pro_1', 'Angiogenic_1', 'Macro_SPP1', 'Ecm_myCAF_1', 'Wound_myCAF', 'Macro_CXCL10']

cellphonedb_clu1 = ['NDRG1+Epi', 'Cycling_Epi', 'Macro_SPP1', 'Macro_CXCL10', 'Ecm_myCAF_1', 'Wound_myCAF', 'Angiogenic_1', 'Angiogenic_2', 'End_progenitor', 'Treg_1']
cellphonedb_clu2 = ['Macro_CXCL10', 'Angiogenic_1', 'Angiogenic_2', 'End_progenitor', 'Cycing_CAF', 'Ecm_myCAF_1', 'Wound_myCAF', 'NDRG1+Epi', 'Cycling_Epi', 'CD8Teff_2', 'Macro_SPP1']

cc = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb_key.csv", index_col=0)
res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_gene_compare.csv", index_col = 0)
grade = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_gene_compare.csv", index_col=0)

for i in cc.index:
    for col in cc.columns:
      if (cc.loc[i, col] > 0):
        idx = i
        flag = False
        if (len(idx.split('_')) > 1):
            flag = True
            idx = idx.split('_')[0]
            gene_c = i.split('_')[1]
        gene_a = idx.split('|')[0]
        gene_b = idx.split('|')[1]
        cluster_1 = col.split('|')[0]
        cluster_2 = col.split('|')[1]
        if ((cluster_1 not in cellphonedb_clu1) | (cluster_2 not in cellphonedb_clu2)):
            continue
        p = pd.DataFrame(data.loc[data.sub_cluster == cluster_1].groupby('sample').mean().loc[:, gene_a])
        p.loc[:, 'type'] = res.loc[p.index, 'type']
        p_t = p.loc[p.type == 'Tumor']
        p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
        p.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/" + gene_a + "-" + cluster_1 + "_type.csv")
        p_t.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/" + gene_a + "-" + cluster_1 + "_grade.csv")
        p = pd.DataFrame(data.loc[data.sub_cluster == cluster_2].groupby('sample').mean().loc[:, gene_b])
        p.loc[:, 'type'] = res.loc[p.index, 'type']
        p_t = p.loc[p.type == 'Tumor']
        p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
        p.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/" + gene_b + "-" + cluster_2 + "_type.csv")
        p_t.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/" + gene_b + "-" + cluster_2 + "_grade.csv")
        if (flag):
            p = pd.DataFrame(data.loc[data.sub_cluster == cluster_2].groupby('sample').mean().loc[:, gene_c])
            p.loc[:, 'type'] = res.loc[p.index, 'type']
            p_t = p.loc[p.type == 'Tumor']
            p_t.loc[:, ['grade', 'stage']] = grade.loc[p_t.index, ['grade', 'stage']]
            p.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/" + gene_c + "-" + cluster_2 + "_type.csv")
            p_t.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/" + gene_c + "-" + cluster_2 + "_grade.csv")


```

# 读取所有的文件
```R
setwd("/public/home/yuwenqi/sc-data/selected/9/module/cci_res/sc_diff/cellphonedb/")

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