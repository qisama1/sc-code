```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
genes = c('SLC1A3','GLUL','SLC1A5','SLC38A1','GLS','GLUD1','OGDH','SUCLG1','GGT6','GGT1','GGT7','GGT5','GSS')
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/gene_df.csv")
```

```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]
major = meta.loc[meta.ident == 'Macro_CXCL10'].index

df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/gene_df.csv", index_col = 0).T
df_use = df.loc[major, ['SLC1A3', 'GLUL', 'SLC1A5', 'SLC38A1', 'GLS', 'GLUD1', 'OGDH', 'SUCLG1', 'GGT7', 'GGT5', 'GSS']]
basic_type = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_compare.csv", index_col = 0)
basic_grade = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_compare.csv", index_col = 0)
meta['type'] = basic_type.loc[meta['sample']].type.values
df_use['type'] = meta.loc[df_use.index, 'type']
df_use.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/CXCL10/gene_plotdata_type.csv")
# meta_t = meta.loc[df_use.loc[df_use.type == 'tumor'].index, ]
# df_use2 = df_use.loc[meta_t.index, ]
# df_use2['grade'] = basic_grade.loc[meta_t['orig.ident']].grade.values
# df_use2['grade'] = df_use2['grade'].replace('grade1&2', 'Low')
# df_use2['grade'] = df_use2['grade'].replace('grade3', 'High')
# df_use2.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/gene_plotdata_grade.csv")

df_use2 = df_use.loc[df_use.type == 'Tumor']
df_use2['grade'] = basic_grade.loc[meta.loc[df_use2.index, 'sample']].stage.values
#df_use2['grade'] = df_use2['grade'].replace('grade1&2', 'Low')
#df_use2['grade'] = df_use2['grade'].replace('grade3', 'High')
df_use2.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/CXCL10/gene_plotdata_grade.csv")

sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/sct.csv", index_col = 0)
scn = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/scn.csv", index_col = 0)
sc = pd.concat([sct, scn])

sc = sc.loc[major, ['M_48']]
#meta['type'] = basic_type.loc[meta['orig.ident']].type.values
sc['type'] = meta.loc[sc.index, 'type']
sc.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/CXCL10/sc_plotdata_type.csv")

sc2 = sc.loc[sc.type == 'Tumor']
sc2['grade'] = basic_grade.loc[meta.loc[sc2.index, 'sample']].stage.values
#sc2['grade'] = sc2['grade'].replace('grade1&2', 'Low')
#sc2['grade'] = sc2['grade'].replace('grade3', 'High')
sc2.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/CXCL10/sc_plotdata_grade.csv")


metabolism_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/metabolism_t.csv", index_col = 0)
metabolism_n = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/metabolism_n.csv", index_col = 0)
metabolism = pd.concat([metabolism_t, metabolism_n])
metabolism = metabolism.loc[major, ['Glutamine']]
#meta['type'] = basic_type.loc[meta['orig.ident']].type.values
metabolism['type'] = meta.loc[metabolism.index, 'type']
metabolism.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/CXCL10/metabolism_plotdata_type.csv")


metabolism2 = metabolism.loc[metabolism.type == 'Tumor']
metabolism2['grade'] = basic_grade.loc[meta.loc[metabolism2.index, 'sample']].stage.values
# metabolism2['grade'] = metabolism2['grade'].replace('grade1&2', 'Low')
# metabolism2['grade'] = metabolism2['grade'].replace('grade3', 'High')
metabolism2.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/esca/CXCL10/metabolism_plotdata_grade.csv")
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
  scale_x_discrete(limits=c('Macro_CXCL10', 'Ap_Epi_2', 'Ap_Epi_1', 'Mes_Epi', 'Mucosal_Epi', 'Cycling_Epi', 'Ap_Epi_3', 'NDRG1+Epi', 'A2M+Epi', 'RPLP1+Epi', 'MAFB+Epi', 'Stress_Epi', 'Mono_3', 'Mono_1', 'Mono_2', 'Macro_Cycing', 'Macro_CXCL10', 'pDC1', 'cDC2', 'mature_cDC', 'Macro_SPP1', 'cDC1', 'Mast',  'Th', 'CD8Teff_1', 'CD8Pro_2', 'Treg_1', 'Treg_3', 'NK_1', 'CD4Tn', 'Treg_2', 'Tfh', 'NK_2', 'CD8Tem_3', 'CD8Tem_1', 'CD8Tem_2', 'CD8Tex', 'CD8Pro_1', 'CD8Teff_2', 'B', 'Fib', 'End')) +
  annotate("text", x = 'Macro_CXCL10', y = 5.5, label='***', alpha=1)
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT2-boxplot/esca/plot.pdf", plot = p1, device = 'pdf', width = 6, height = 4, units = 'in')



```