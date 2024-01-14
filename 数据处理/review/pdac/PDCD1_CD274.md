# 读数
```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", row.names=1)

gene = c('PDCD1', 'CD274')
write.csv(as.matrix(scRNA[gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/gene/pdac_gene_df.csv")
```

# vinlin分析
```py
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/gene/pdac_gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col = 0)
meta.loc[meta['sample'].str.contains("N"), 'Tissue'] = 'Normal'
meta.loc[~meta['sample'].str.contains("N"), 'Tissue'] = 'Tumor'


gene_df = gene_df.loc[meta.index,:]
gene_df['cluster'] = meta['cluster']

gene_df.to_csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/violin/pdac_data.csv")
```
# violinplot
```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/violin/pdac_data.csv")
p1 = ggplot(data,aes(x=cluster,y=PDCD1),fill=cluster)+
    scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                '#D9D9D9', '#EEAEEE'))+
    geom_violin(aes(fill=cluster),show.legend = T)+ #outlier.shape = NA
    theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(),
            axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
            plot.title = element_text(hjust=0.5),)+
    labs(x=NULL,
        y=NULL,
        title = "PDCD1")+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=cluster), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['PDCD1']) - 0.15)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
    
                        
p1=p1+theme(legend.position = "none",legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.5, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 14) )
ggsave("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/violin/pdac_data_PDCD1.pdf", width = 2, height = 2.6, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/violin/pdac_data.csv")
p1 = ggplot(data,aes(x=cluster,y=CD274),fill=cluster)+
    scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                '#D9D9D9', '#EEAEEE'))+
    geom_violin(aes(fill=cluster),show.legend = T)+ #outlier.shape = NA
    theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(),
            axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
            plot.title = element_text(hjust=0.5),)+
    labs(x=NULL,
        y=NULL,
        title = "CD274")+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=cluster), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['CD274']) - 0.15)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.position = "none",legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.5, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 14) )
ggsave("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/violin/pdac_data_CD274.pdf", width =2, height = 2.6, units = 'in')

```

# tumor_normal分析
```py
import pandas as pd
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/gene/pdac_gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col = 0)
meta.loc[meta['sample'].str.contains("N"), 'Tissue'] = 'Normal'
meta.loc[~meta['sample'].str.contains("N"), 'Tissue'] = 'Tumor'


gene_df = gene_df.loc[meta.index,:]
clusters = ['Mye', 'TNK', 'End', 'Fib', 'Epi']
for cluster in clusters:
    t = gene_df
    t['type'] = meta['Tissue']
    t = t.loc[meta[meta.cluster == cluster].index]
    t.to_csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/boxplot/" + cluster + "_pdac.csv")   
```

# boxplot_tumor_normal
```R
setwd("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/boxplot/")
fs = list.files('./', '_pdac.csv')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    data = read.csv(x, row.names=1)
    data = melt(data, id.vars = c('type'))

    p1 = ggplot(data,aes(x=variable,y=value),fill=type)+
        scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                    '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                    '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                    '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                    '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                    '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                    '#D9D9D9', '#EEAEEE'))+
        geom_boxplot(aes(fill=type),show.legend = T)+ #outlier.shape = NA
        theme(panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(),
                axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
                plot.title = element_text(hjust=0.5),)+
        labs(x=NULL,
            y=NULL,
            title = NULL)+
        #theme_minimal() 
        #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
        #coord_cartesian(ylim = c(0, 2))+
        stat_compare_means(aes(group=type), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                            label = "p.signif", method='wilcox.test', label.x = 0.55, label.y = max(data['value']) - 0.15)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                            
    p1=p1+theme(legend.position = "none",legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.5, "cm"), legend.title = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 14) )
    
    ggsave(paste0("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/boxplot/plot/",filename, ".pdf"), width = 3, height = 2.6, units = 'in')

})
```

# sumTME_process
```py
import pandas as pd

gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/gene/pdac_gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col = 0)
meta.loc[meta['sample'].str.contains("N"), 'Tissue'] = 'Normal'
meta.loc[~meta['sample'].str.contains("N"), 'Tissue'] = 'Tumor'



gene_df = gene_df.loc[meta.index, ]
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/module.csv")



res = pd.DataFrame()

for m in ['module1', 'module2', 'module3', 'module4']:
    for cluster in module[m].dropna():
        res = pd.concat([res, pd.DataFrame([[cluster, m]], index = [cluster], columns = ['sub_cluster', 'subTME'])])

for pathway in gene_df.columns:
    for idx in res.index:
        cur = meta.loc[meta.ident == idx]
        res.loc[idx, pathway] = gene_df.loc[cur.index, pathway].mean()
res = res.iloc[:, 1:]

res = res.replace('module1', 'subTME-PSE')
res = res.replace('module2', 'subTME-PSE')
res = res.replace('module3', 'subTME-IS+MRM')
res = res.replace('module4', 'subTME-MRM')

res.to_csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/subTME_Diff/pdac_plotdata.csv")
```

# sumTME_plot

```R
library(ggsci)
library(ggplot2)
library(cowplot)
library(data.table)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(tidyverse)

data = read.csv("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/subTME_Diff/pdac_plotdata.csv", row.names=1)
data = melt(data, id.vars = c('subTME'))

p1 = ggplot(data,aes(x=variable,y=value),fill=subTME)+
    scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                '#D9D9D9', '#EEAEEE'))+
    geom_boxplot(aes(fill=subTME),show.legend = T)+ #outlier.shape = NA
    theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(),
            axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
            plot.title = element_text(hjust=0.5),)+
    labs(x=NULL,
        y=NULL,
        title = NULL)+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='wilcox.test', label.x = 0.55, label.y = max(data['value']) - 0.1,size=5)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.position = "none",legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.5, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 20) )

ggsave("/public/home/yuwenqi/sc-data/selected/review/PDCD1_CD274/subTME_Diff/plot/pdac_plot.pdf", width = 2, height = 2.6, units = 'in')


```