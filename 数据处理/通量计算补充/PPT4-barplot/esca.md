```py
barplot_use = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/barplot_meta.csv")
metabolism = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/esca/metabolism_t.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]

metabolism['ident'] = meta_t.loc[metabolism.index, 'ident']

# Epi
major = meta_t.loc[meta_t.cluster == 'Epi'].index
metabolism.loc[major, barplot_use.Epi.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/esca/esca_barplot_epi.csv")
# Mye
major = meta_t.loc[meta_t.cluster == 'Mye'].index
metabolism.loc[major, barplot_use.Mye.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/esca/esca_barplot_Mye.csv")
# APOE
major = meta_t.loc[meta_t.ident == 'Macro_SPP1_APOE'].index
metabolism.loc[major, barplot_use.APOE.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/esca/esca_barplot_Macro_SPP1_APOE.csv")
# CXCL10
major = meta_t.loc[meta_t.ident == 'Macro_CXCL10'].index
metabolism.loc[major, barplot_use.APOE.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/esca/esca_barplot_Macro_CXCL10.csv")
```