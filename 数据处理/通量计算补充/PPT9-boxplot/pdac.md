```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
genes = c('SLC1A3', 'GLUL', 'SLC1A5', 'SLC38A1', 'GLS', 'GLUD1', 'OGDH' , 'SUCLG1', 'GGT6', 'GGT1', 'GGT7', 'GGT5')
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT9-boxplot/pdac/gene_df.csv")
```

```python
import pandas as pd
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]
major1 = meta_t.loc[meta_t.ident == 'Macro_APOE'].index
major2 = meta_t.loc[meta_t.cluster == 'Epi'].index

df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT9-boxplot/pdac/gene_df.csv", index_col = 0).T
df_use1 = df.loc[major1, ['SLC1A3', 'GLUL', 'SLC1A5', 'SLC38A1', 'GLS', 'GLUD1', 'OGDH' , 'SUCLG1', 'GGT6', 'GGT1', 'GGT7', 'GGT5']]
df_use2 = df.loc[major2, ['SLC1A3', 'GLUL', 'SLC1A5', 'SLC38A1', 'GLS', 'GLUD1', 'OGDH' , 'SUCLG1', 'GGT6', 'GGT1', 'GGT7', 'GGT5']]
df_use1['type'] = 'Macro_SPP1_APOE'
df_use2['type'] = 'Tumor cell'
pd.concat([df_use1, df_use2]).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT9-boxplot/pdac/plot_data.csv")
```

```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT2-boxplot/pdac/plot_data.csv", row.names = 1)
p1 = ggplot(data,aes(x=ident,y=GLUL),fill=ident)+
  scale_fill_manual(values=c('#FFDEAD','#7CFC00','#20B2AA',  '#F0F8FF', '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', '#C2C2C2', '#C1CDC1', '#CD1076', '#B03060', '#AB82FF', '#8E8E38', '#FF1493', '#708090', '#4F94CD', '#473C8B', '#212121', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#FF1493', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95', '#D9D9D9', '#EEAEEE', '#D9D9D9', '#EEAEEE', '#FF1493', '#708090', '#4F94CD', '#473C8B', '#212121', '#00E5EE', '#006400', '#8B0000', '#8B5F65'))+
  geom_boxplot(aes(fill=ident),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL, 
       title = "GLUL") +
  #annotate("text", label='***', x = 'Macro_APOE/CTSZ', y = 4)
  #theme_minimal() 
  #geom_sina(aes(color=cluster), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 0.5))+
  scale_x_discrete(limits=c('Macro_APOE', 'Macro_SPP1_APOE',  'Basal', 'Ductal_cell_A_5', 'Mesenchymal', 'Endocrine_cell', 'Ductal_cell_A_3', 'Ductal_cell_B_4', 'Interferon_Epi', 'Ductal_cell_A_4', 'Ductal_cell_B_1', 'Acinar_cell_1', 'Ductal_cell_A_2', 'Ductal_cell_B_5', 'Acinar_cell_3', 'Ductal_cell_B_3', 'Ductal_cell_A_1', 'Ductal_cell_B_2', 'Acinar_cell_2', 'TRM', 'Mono_1', 'Mono_2',  'Macro_Cycing', 'Macro_CXCL10', 'cDC2', 'pDC1', 'Macro_ACP5', 'mature_cDC', 'Macro_SPP1', 'Macro_CREM', 'cDC1', 'Mast', 'Langerhans_like',  'Th', 'CD4Tn', 'CD8Pro', 'CD4Tex', 'Treg', 'CD8Tem_1', 'CD8Tem_2', 'CD8Tex', 'CD4Tem_1', 'T_progenitors', 'B', 'Fib', 'End')) +
  annotate("text", x = 'Macro_APOE', y = 4.5, label='***', alpha=1)
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT2-boxplot/pdac/plot.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')

```