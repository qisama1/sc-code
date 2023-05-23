```py
barplot_use = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/barplot_meta.csv")
metabolism = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/pdac/metabolism_t.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]

metabolism['ident'] = meta_t.loc[metabolism.index, 'ident']

# Epi
major = meta_t.loc[meta_t.cluster == 'Epi'].index
metabolism.loc[major, barplot_use.Epi.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/pdac/pdac_barplot_epi.csv")
# Mye
major = meta_t.loc[meta_t.cluster == 'Mye'].index
metabolism.loc[major, barplot_use.Mye.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/pdac/pdac_barplot_mye.csv")
# APOE
major = meta_t.loc[meta_t.ident == 'Macro_APOE'].index
metabolism.loc[major, barplot_use.APOE.dropna()].mean().to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/pdac/pdac_barplot_Macro_APOE.csv")
```