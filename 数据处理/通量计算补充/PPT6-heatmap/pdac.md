```py
data_use = pd.read_csv("/public/home/yuwenqi/data/Data26/ppt13-heatmap/part.csv")
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT6-heatmap/module.csv", index_col = 0)
sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/pdac/sct.csv", index_col = 0)
gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/all2.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]

major = meta_t.loc[meta_t.ident.str.contains('Macro')].index

t = sct.loc[major, data_use.B.dropna()]
t.columns = module.loc[t.columns, 'name']
gsva = gsva.loc[gsva.index.isin(major), ['M1 Macrophage Polarization', 'Pro-inflammatory in myeloid cells', 'M2 Macrophage Polarization', 'Anti-inflammatory in myeloid cells', 'Angiogenesis']]
gsva['ident'] = meta.loc[gsva.index, 'ident']

res = pd.DataFrame(index = t.columns, columns = gsva.columns)
for i in res.index:
    for j in res.columns:
        cur = spearmanr(t[i], gsva[j])
        if (cur[1] > 0.05):
            res.loc[i, j] = 0
        else :
            res.loc[i, j] = cur[0]

res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT6-heatmap/pdac/plot_data.csv")
```