# 按照差异基因选
```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/mye_info.csv", index_col = 0)
mye_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/mye/mye_res.csv", index_col=0)
marker = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/mye/myeloid-res0.8.csv")

sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.ident))
sg = sg.fillna(0)

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, mye_res.loc[mye_res.cluster == t.cluster].index] = score

for i in sg.columns:
    t = sg.loc[:, i]
    pd.DataFrame(t.sort_values(ascending = False)).to_csv("/public/home/yuwenqi/sc-data/selected/35/mye/35_mye_score/"+ i + "_pct.csv")
```