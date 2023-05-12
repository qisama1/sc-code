```python
import pandas as pd
tcga = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/tcga_tpm_log.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ESCA_sample_meta2.csv", index_col = 0).T
gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/needed_gene.csv")
tcga_use = tcga.loc[meta.index, ]
res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_type.csv", index_col=0)
tcga_use.loc[:, 'type'] = res.loc[tcga_use.index, 'type']     
os_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/os/os_meta.csv", index_col = 0)
res_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_stage.csv", index_col=0)
res_t2 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_grade.csv", index_col=0)


gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ana_data.csv")
gene_sub = gene.loc[gene.type == 'ESCA']
for idx in gene_sub.index:
    v = gene_sub.loc[idx, ]
    gene_a = v.gene.split('|')[0]
    gene_b = v.gene.split('|')[1]
    cluster_1 = v.cluster.split('|')[0]
    cluster_2 = v.cluster.split('|')[1]
    p = pd.DataFrame(tcga_use.loc[:, [gene_a]])
    p2 = pd.DataFrame(tcga_use.loc[:, [gene_a]])
    p2 = p2[~p2.index.duplicated()]
    p2.index = p2.index.str[0:12]
    p_t = p2.loc[res_t.index,]
    p_t2 = p2.loc[res_t2.index,]
    p.loc[:, 'type'] = res.loc[p.index, 'type']
    p_t.loc[:, ['stage', 'stage2']] = res_t.loc[p_t.index, ['stage', 'stage2']]
    p_t2.loc[:, ['grade', 'grade2']] = res_t2.loc[p_t2.index, ['grade', 'grade2']]
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/" + gene_a + "_type.csv")
    p_t.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/" + gene_a + "_stage.csv")
    p_t2.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/" + gene_a + "_grade.csv")
    p = pd.DataFrame(tcga_use.loc[:, [gene_b]])
    p2 = pd.DataFrame(tcga_use.loc[:, [gene_b]])
    p2 = p2[~p2.index.duplicated()]
    p2.index = p2.index.str[0:12]
    p_t = p2.loc[res_t.index,]
    p_t2 = p2.loc[res_t2.index,]
    p.loc[:, 'type'] = res.loc[p.index, 'type']
    p_t.loc[:, ['stage', 'stage2']] = res_t.loc[p_t.index, ['stage', 'stage2']]
    p_t2.loc[:, ['grade', 'grade2']] = res_t2.loc[p_t2.index, ['grade', 'grade2']]
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/" + gene_b + "_type.csv")
    p_t.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/" + gene_b + "_stage.csv")
    p_t2.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/" + gene_b + "_grade.csv")
        
gene_os = pd.read_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/gene_os.csv")
for gene in gene_os.gene:
    p = pd.DataFrame(tcga_use.loc[:, [gene]])
    p = pd.concat([p, os_res], axis=1)
    p.to_csv("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/os/" + gene + "_os.csv")

```

```R
setwd("/public/home/yuwenqi/sc-data/selected/target-analysis-0403/ESCA/tcga_diff/")
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

fs = list.files('./', 'stage')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x)
    data=melt(data, id.vars=c('stage', 'stage2', 'X'))
    boxplot_type3(data, paste0("./plot3/", filename, "_stage2.pdf"))
})

lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x)
    data=melt(data, id.vars=c('stage', 'stage2', 'X'))
    boxplot_type2(data, paste0("./plot2/", filename, "_stage.pdf"))
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



fs = list.files('./', 'grade')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x)
    data=melt(data, id.vars=c('grade', 'grade2', 'X'))
    boxplot_type5(data, paste0("./plot5/", filename, "_grade2.pdf"))
})

lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x)
    data=melt(data, id.vars=c('grade', 'grade2', 'X'))
    boxplot_type4(data, paste0("./plot4/", filename, "_grade.pdf"))
})

boxplot_type4 = function(data, filename) {
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
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=grade),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}

boxplot_type5 = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"grade2"])))]
    p=ggboxplot(data, x="variable", y="value", color = "grade2",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=grade2),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=grade2),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}
```