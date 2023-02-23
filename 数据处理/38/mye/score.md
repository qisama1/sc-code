# 取基因表达
```R
scRNA = qread("...")
write.csv(as.matrix(scRNA@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/38/score/mye/mye_df.csv")
```

# 计算单个簇的基因分数

```python
df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/score/mye/mye_df.csv", index_col=0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/mye_info2.csv", index_col = 0)
mye_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/score/mye/mye_res.csv", index_col=0)

df = df.T
cnt = pd.DataFrame(index = df.columns, columns = set(meta.ident))
cnt = cnt.fillna(0)

for i in range(len(df)):
    t = df.iloc[i,]
    score = df.loc[meta.loc[meta.ident == t.name].index].sum() / 4
    cnt.loc[t.loc[t > score].index, meta.loc[t.name, 'ident']] += 1

sg = pd.DataFrame(index = cnt.index, columns = cnt.columns)

for i in cnt.columns:
    t = cnt.loc[:, i]
    for j in cnt.index:
        if (cnt.loc[j, i] == 0):
            sg.loc[j, i] = 0
            continue
        p_cur = cnt.loc[j, i] / len(meta.loc[meta.ident == i])
        p_other = (cnt.loc[j, ].sum() - cnt.loc[j, i]) / len(meta.loc[meta.ident != i])
        sg.loc[j, i] = 1 - (float(p_other) / float(p_cur))

marker = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/mye/myeloid-res0.8.csv")

# 去除RPL、RPS、MT- 开头的这些基因
sg = sg.loc[(sg.index.str[0:3] != 'RPL') & (sg.index.str[0:3] != 'RPS') & (sg.index.str[0:2] != 'MT')]
marker = marker.loc[(marker.gene.str[0:3] != 'RPL') & (marker.gene.str[0:3] != 'RPS') & (marker.gene.str[0:2] != 'MT')]
for i in sg.columns:
    t = sg.loc[:, i]
    cur = pd.DataFrame(sg.index[t > 0.4], columns = ['gene'])
    mark = set(marker.loc[marker.cluster == mye_res.loc[i,'cluster']].loc[marker.p_val <= 0.05].loc[marker.avg_log2FC >= 0.5].gene)
    pd.DataFrame(set(cur['gene']) & mark).to_csv("/public/home/yuwenqi/sc-data/selected/38/score/sg" + i + "_0.4.csv")

# 大于其他的最大值
for i in sg.columns:
    t = sg.loc[:, i]
    cur = pd.DataFrame(sg.index[t > sg.drop(i,axis=1).max(axis=1)], columns = ['gene'])
    mark = set(marker.loc[marker.cluster == mye_res.loc[i,'cluster']].loc[marker.p_val <= 0.05].loc[marker.avg_log2FC >= 0.5].gene)
    pd.DataFrame(set(cur['gene']) & mark).to_csv("/public/home/yuwenqi/sc-data/selected/38/score/othermax/sg" + i + "_othermax.csv")

# 大于0.1 选top20，并去除RPL、RPS、MT- 开头的这些基因
for i in sg.columns:
    if (type(mye_res.loc[i, 'cluster']) == np.int64) :
        mark = set(marker.loc[marker.cluster == (mye_res.loc[i, 'cluster'])].loc[marker.p_val <= 0.05].loc[marker.avg_log2FC >= 0.5].gene)
    else :
        mark = set(marker.loc[marker.cluster.isin (mye_res.loc[i, 'cluster'].values)].loc[marker.p_val <= 0.05].loc[marker.avg_log2FC >= 0.5].gene)
    t = sg.loc[:, i]
    t = t.loc[mark,]
    t = t.loc[t > 0.1]
    cur = pd.DataFrame(t.sort_values()[-20:].index, columns = ['gene'])
    cur.to_csv("/public/home/yuwenqi/sc-data/selected/38/score/mye/sg" + i + "_0.1_top20.csv")

# 结合marker
CD4_Tn = pd.read_csv("/public/home/yuwenqi/sc-data/selected/38/score/sgCD4+Tn.csv", index_col=0)
marker_gene = set(marker.loc[marker.cluster == 0].loc[marker.p_val <= 0.05].loc[marker.avg_log2FC >= 0.5].gene)
set(CD4_Tn['0']) & marker_gene

```


# 按照差异基因选
```python
sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.ident))
sg = sg.fillna(0)

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, mye_res.loc[mye_res.cluster == t.cluster].index] = score

for i in sg.columns:
    t = sg.loc[:, i]
    pd.DataFrame(t.sort_values(ascending = False)).to_csv("/public/home/yuwenqi/sc-data/selected/38/score/mye/diff_mye/sg" + i + "_pct.csv")
```