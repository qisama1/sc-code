# 按照差异基因选
```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/epi_info.csv", index_col = 0)
epi_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/score/epi/epi_res.csv", index_col=0)
marker = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/epi/epi-res0.3-harmony.csv")

sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.ident))
sg = sg.fillna(0)
sg = sg.loc[(sg.index.str[0:3] != 'RPL') & (sg.index.str[0:3] != 'RPS') & (sg.index.str[0:3] != 'MT-') & (sg.index.str[0:3] != 'RLS')]

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, epi_res.loc[epi_res.cluster == t.cluster].index] = score

for i in sg.columns:
    t = sg.loc[:, i]
    res = pd.DataFrame(t.loc[t > 0.3].sort_values(ascending = False)[0:30])
    res.to_csv("/public/home/yuwenqi/sc-data/selected/38/score/epi/diff_epi/sg" + i + "_pct.csv")
```