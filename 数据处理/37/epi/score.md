# 按照差异基因选
```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/epi_info.csv", index_col = 0)
epi_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/epi/epi_res.csv", index_col=0)
marker = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/epi/epi-res0.5-harmony.csv")

sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.ident))
sg = sg.fillna(0)

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, epi_res.loc[epi_res.cluster == t.cluster].index] = score

for i in sg.columns:
    t = sg.loc[:, i]
    pd.DataFrame(t.sort_values(ascending = False)).to_csv("/public/home/yuwenqi/sc-data/selected/37/epi/37_epi_score/"+ i + "_pct.csv")
```