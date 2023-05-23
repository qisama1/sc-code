```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
genes = c('SLC1A3','GLUL','SLC1A5','SLC38A1','GLS','GLUD1','OGDH','SUCLG1','GGT6','GGT1','GGT7','GGT5','GSS')
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/bc/gene_df.csv")
```

```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
meta_t = meta.loc[meta.Tissue == 'Tumor']
df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/bc/SLC1A3_df.csv", index_col = 0).T
df = df.loc[meta_t.index, ]
df['ident'] = meta_t['sub_cluster']

# get Epi、Mye、TNK
df.loc[meta_t[~meta_t.cluster.isin(['Epi', 'Mye', 'TNK'])].index, 'ident'] = meta_t.cluster
df.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/bc/SLC1A3_diff.csv")

sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/bc/sct2.csv", index_col = 0)
sct = sct.loc[meta_t.index, ['M_48', 'M_25', 'M_26']]
sct['ident'] = meta_t['sub_cluster']
sct.loc[meta_t[~meta_t.cluster.isin(['Epi', 'Mye', 'TNK'])].index, 'ident'] = meta_t.cluster
sct.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/bc/module_diff.csv")


metabolism = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/bc/metabolism_t2.csv", index_col = 0)
metabolism = metabolism.loc[meta_t.index, ['Glutamate', 'Glutamine']]
metabolism['ident'] = meta_t['sub_cluster']
metabolism.loc[meta_t[~meta_t.cluster.isin(['Epi', 'Mye', 'TNK'])].index, 'ident'] = meta_t.cluster
metabolism.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/bc/metabolism_diff.csv")


set(meta_t.loc[meta_t.cluster == 'Epi'].ident)
set(meta_t.loc[meta_t.cluster == 'Mye'].ident)
set(meta_t.loc[meta_t.cluster == 'TNK'].ident)
set(meta_t.cluster)
```

```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/bc/SLC1A3_diff.csv", row.names = 1)
p1 = ggplot(data,aes(x=ident,y=SLC1A3),fill=ident)+
  scale_fill_manual(values=c('#FFDEAD','#7CFC00','#20B2AA',  '#F0F8FF', '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', '#C2C2C2', '#C1CDC1', '#CD1076', '#B03060', '#AB82FF', '#8E8E38', '#FF1493', '#708090', '#4F94CD', '#473C8B', '#212121', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95', '#D9D9D9', '#EEAEEE', '#FF1493', '#708090', '#4F94CD', '#473C8B', '#212121', '#00E5EE', '#006400', '#8B0000', '#8B5F65', '#CD8C95', '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=ident),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL, 
       title = "SLC1A3") +
  #annotate("text", label='***', x = 'Macro_APOE/CTSZ', y = 4)
  #theme_minimal() 
  #geom_sina(aes(color=cluster), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 0.5))+
  scale_x_discrete(limits=c('Macro_APOE', 'Basal', 'Cycling', 'ML_5', 'LP_1', 'ML_3', 'ML_2', 'ML_4', 'Her2+', 'ML_6', 'ML_1','Mono_1', 'Mono_2', 'DC_2', 'DC_1', 'Macro_CX3CR1', 'Macro_SPP1', 'Intermediate', 'NK', 'T_helper_1', 'T_helper_2', 'CD4Tn', 'Treg', 'CD8Teff', 'gd_T', 'CD8Tem_1', 'CD8Tem_2', 'ProT', 'End', 'Fib')) +
  annotate("text", x = 'Macro_APOE', y =3.3, label='***', alpha=1)
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT7-boxplot/bc/SLC1A3_diff.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')
```