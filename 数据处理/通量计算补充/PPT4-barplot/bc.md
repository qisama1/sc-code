```py
barplot_use = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/barplot_meta.csv")
metabolism = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/bc/metabolism_t2.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
meta_t = meta.loc[meta.Tissue == 'Tumor']

metabolism.index = metabolism.index.str.replace('.', '-')
meta_t['ident'] = meta_t['sub_cluster']
metabolism = metabolism.loc[meta_t.index, ]
metabolism['ident'] = meta_t.loc[metabolism.index, 'ident']

# Epi
major = meta_t.loc[meta_t.cluster == 'Epi'].index
metabolism.loc[major, barplot_use.Epi.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/bc/bc_barplot_epi.csv")
# Mye
major = meta_t.loc[meta_t.cluster == 'Mye'].index
metabolism.loc[major, barplot_use.Mye.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/bc/bc_barplot_Mye.csv")
# APOE
major = meta_t.loc[meta_t.ident == 'Macro_APOE'].index
metabolism.loc[major, barplot_use.APOE.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/bc/bc_barplot_Macro_APOE.csv")

```