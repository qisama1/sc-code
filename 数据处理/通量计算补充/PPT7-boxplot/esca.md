```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
write.csv(as.matrix(scRNA['SLC1A3', ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/esca/SLC1A3_df.csv")
```

```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]
df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/esca/SLC1A3_df.csv", index_col = 0).T
df = df.loc[meta_t.index, ]
df['ident'] = meta_t['ident']

# get Epi、Mye、TNK
df.loc[meta_t[~meta_t.cluster.isin(['Epi', 'Mye', 'TNK'])].index, 'ident'] = meta_t.cluster
df.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/esca/SLC1A3_diff.csv")

sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/sct.csv", index_col = 0)
sct = sct.loc[meta_t.index, ['M_48', 'M_25', 'M_26']]
sct['ident'] = meta_t['ident']
sct.loc[meta_t[~meta_t.cluster.isin(['Epi', 'Mye', 'TNK'])].index, 'ident'] = meta_t.cluster
sct.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/esca/module_diff.csv")


metabolism = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/metabolism_t.csv", index_col = 0)
metabolism = metabolism.loc[meta_t.index, ['Glutamate', 'Glutamine']]
metabolism['ident'] = meta_t['ident']
metabolism.loc[meta_t[~meta_t.cluster.isin(['Epi', 'Mye', 'TNK'])].index, 'ident'] = meta_t.cluster
metabolism.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/esca/metabolism_diff.csv")

set(meta_t.loc[meta_t.cluster == 'Epi'].ident)
set(meta_t.loc[meta_t.cluster == 'Mye'].ident)
set(meta_t.loc[meta_t.cluster == 'TNK'].ident)
set(meta_t.cluster)
```

```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT2-boxplot/esca/plot_data.csv", row.names = 1)
p1 = ggplot(data,aes(x=ident,y=GLUL),fill=ident)+
  scale_fill_manual(values=c('#FFDEAD','#7CFC00','#20B2AA',  '#F0F8FF', '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', '#C2C2C2', '#C1CDC1', '#CD1076', '#B03060', '#AB82FF', '#8E8E38', '#FF1493', '#708090', '#4F94CD', '#473C8B', '#212121', '#00E5EE', '#006400', '#8B0000', '#FF1493', '#CD8C95', '#D9D9D9', '#EEAEEE', '#8B5F65', '#708090', '#4F94CD', '#473C8B', '#8B5F65', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95', '#D9D9D9', '#EEAEEE', '#D9D9D9', '#EEAEEE', '#FF1493', '#708090', '#4F94CD', '#473C8B', '#212121', '#00E5EE', '#006400', '#8B0000', '#8B5F65'))+
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
  scale_x_discrete(limits=c('Macro_SPP1_APOE', 'Ap_Epi_2', 'Ap_Epi_1', 'Mes_Epi', 'Mucosal_Epi', 'Cycling_Epi', 'Ap_Epi_3', 'NDRG1+Epi', 'A2M+Epi', 'RPLP1+Epi', 'MAFB+Epi', 'Stress_Epi', 'Mono_3', 'Mono_1', 'Mono_2', 'Macro_Cycing', 'Macro_CXCL10', 'pDC1', 'cDC2', 'mature_cDC', 'Macro_SPP1', 'cDC1', 'Mast',  'Th', 'CD8Teff_1', 'CD8Pro_2', 'Treg_1', 'Treg_3', 'NK_1', 'CD4Tn', 'Treg_2', 'Tfh', 'NK_2', 'CD8Tem_3', 'CD8Tem_1', 'CD8Tem_2', 'CD8Tex', 'CD8Pro_1', 'CD8Teff_2', 'B', 'Fib', 'End')) +
  annotate("text", x = 'Macro_SPP1_APOE', y = 5.5, label='***', alpha=1)
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT2-boxplot/esca/plot.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')



```