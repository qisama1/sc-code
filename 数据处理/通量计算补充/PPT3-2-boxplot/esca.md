```R
genes = c('IL1B', 'CD40', 'IL6', 'CD86', 'APOE', 'CD163', 'CSF1R', 'CTSZ', 'GPNMB', 'CCL18')
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/gene_df.csv")
```


```py
df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]
meta_use = meta_t.loc[meta_t.cluster == 'Mye']
df_use = df.loc[meta_use.index, ]
df_use['ident'] = meta_t['ident']
df_use.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv")
```

```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))

bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="IL1B", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "IL1B")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +

  annotate("text", x = 'Macro_SPP1_APOE', y =4.8, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/IL1B.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="CD40", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CD40")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  annotate("text", x = 'Macro_SPP1_APOE', y =4.8, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/CD40.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="IL6", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "IL6")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 2.5, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/IL6.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="CD86", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CD86")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 2.5, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/CD86.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="APOE", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "APOE")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 7, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/APOE.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="CD163", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CD163")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 4, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/CD163.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="CSF1R", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CSF1R")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 4, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/CSF1R.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="CTSZ", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CTSZ")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 4, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/CTSZ.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="GPNMB", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "GPNMB")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 4, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/GPNMB.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/plot_data.csv", row.names = 1)
data$ident <- factor(data$ident,levels = c('Macro_SPP1_APOE', 'cDC2', 'Macro_Cycing', 'Mono_1', 'Mono_3', 'Macro_SPP1', 'pDC1', 'cDC1', 'Mono_2', 'Macro_CXCL10', 'Mast', 'mature_cDC'))
bioCol=c("#ff3399", "#1c79c0",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95')
bioCol=bioCol[1:length(levels(factor(data[,"ident"])))]
p=ggboxplot(data, x="ident", y="CCL18", color = "ident",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CCL18")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=ident),position = position_jitterdodge(),size=1.0) +
  
  annotate("text", x = 'Macro_SPP1_APOE', y = 4, label='***', alpha=1) +
  theme(legend.position = 'none')

ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT3-2-boxplot/esca/CCL18.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')
```