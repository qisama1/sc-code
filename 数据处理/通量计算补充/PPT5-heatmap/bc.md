```py
module_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT5-heatmap/modules_epi.csv", index_col = 0)
module_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT5-heatmap/modules_mye.csv", index_col = 0)

sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/bc/sct.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
meta_t = meta.loc[meta.Tissue == 'Tumor']
meta_t['ident'] = meta_t['sub_cluster']
# Epi
major = meta_t.loc[meta_t.cluster == 'Epi'].index
t = sct.loc[major, module_epi.A]
t['ident'] = meta_t.loc[t.index, 'ident']
t.groupby('ident').mean().T.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT5-heatmap/bc/bc_heatmap_epi.csv")
# Mye
major = meta_t.loc[meta.cluster == 'Mye'].index
t = sct.loc[major, module_mye.A]
t['ident'] = meta_t.loc[t.index, 'ident']
t.groupby('ident').mean().T.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT5-heatmap/bc/bc_heatmap_mye.csv")
```