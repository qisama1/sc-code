```python
import pandas as pd
tcga = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/tcga_tpm_log.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/BRCA_sample_meta2.csv", index_col = 0).T
tcga.index = tcga.index.str[2:-3]
gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/TCGA/tcga_part1_gene.csv")
tcga_use = tcga.loc[meta.index, gene.Gene]
res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/module/tpm_compare/tpm_module_type.csv", index_col=0)
tcga_use.loc[:, 'type'] = res.loc[tcga_use.index, 'type']     
os_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/os/os_res.csv", index_col = 0)
res_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/module/tpm_compare/tpm_module_stage.csv", index_col=0)

cc = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cellphonedb-key.csv", index_col = 0)
for idx in cc.index:
    for col in cc.columns:
      if (cc.loc[idx, col] > 0):
        gene_a = idx.split('|')[0]
        gene_b = idx.split('|')[1]
        if (gene_b == 'ITGA5_ITGB1'):
            continue
        p = pd.DataFrame(tcga_use.loc[:, [gene_a]])
        p_t = p.loc[res_t.index,]
        p.loc[:, 'type'] = res.loc[p.index, 'type']
        p_t.loc[:, ['stage', 'stage2']] = res_t.loc[p_t.index, ['stage', 'stage2']]
        p.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/diff/cellphonedb/" + gene_a + "_type.csv")
        p_t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/diff/cellphonedb/" + gene_a + "_grade.csv")
        p = pd.DataFrame(tcga_use.loc[:, [gene_b]])
        p_t = p.loc[res_t.index,]
        p.loc[:, 'type'] = res.loc[p.index, 'type']
        p_t.loc[:, ['stage', 'stage2']] = res_t.loc[p_t.index, ['stage', 'stage2']]
        p.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/diff/cellphonedb/" + gene_b + "_type.csv")
        p_t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/diff/cellphonedb/" + gene_b + "_grade.csv")

for idx in cc.index:
    for col in cc.columns:
      if (cc.loc[idx, col] > 0):
        gene_a = idx.split('|')[0]
        gene_b = idx.split('|')[1]
        if (gene_b == 'ITGA5_ITGB1'):
            continue
        p = pd.DataFrame(tcga_use.loc[:, [gene_a]])
        p = pd.concat([p, os_res], axis=1)
        p.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/os/" + gene_a + "_os.csv")
        p = pd.DataFrame(tcga_use.loc[:, [gene_b]])
        p = pd.concat([p, os_res], axis=1)
        p.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/os/" + gene_b + "_os.csv")
```

```R
setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/result/tcga_res/diff/cellphonedb/")
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
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('stage', 'stage2'))
    boxplot_type3(data, paste0("./plot3/", filename, "_stage2.pdf"))
})

lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('stage', 'stage2'))
    boxplot_type3(data, paste0("./plot2/", filename, "_stage.pdf"))
})

boxplot_type2 = function(data, filename) {
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
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=stage),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}

boxplot_type3 = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"stage2"])))]
    p=ggboxplot(data, x="variable", y="value", color = "stage2",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=stage2),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=stage2),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}
```