# Mye重新改造
```R
Idents(mye) = mye$seurat_clusters
new.cluster.ids <- c("Macro_APOE", "Mono_1", "Mono_2", "Macro_FTL", "DC_1", "Macro_Pro", "Macro_CXCL8", "Macro_CCL19", "DC_2")

names(new.cluster.ids) <- levels(mye)
mye = RenameIdents(mye, new.cluster.ids)
mye[['sub_cluster']] = Idents(mye)
mye = subset(mye, ident != 'del')

cellinfo <- subset(mye@meta.data, select = c("orig.ident", "sub_cluster", 'dataset'))
write.csv(cellinfo, "/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/mye_meta2.csv")
```
# 整合meta
```python
mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/mye_meta2.csv", index_col=0)
tnk = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/t_meta.csv", index_col=0)
tnk = tnk.drop('MMRStatus', axis=1)
end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/End/end_meta.csv", index_col=0)
fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Fib/fib_meta2.csv", index_col=0)
epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Epi2/epi_meta.csv", index_col=0)
epi['sub_cluster'] = epi['ident']
```

# 计算比例
```python
mye['cluster'] = 'Mye'
tnk['cluster'] = 'Mye'
epi['cluster'] = 'Epi'
fib['cluster'] = 'Fib'
end['cluster'] = 'End'

meta = pd.concat([mye, tnk, epi, fib, end])
meta = meta.loc[meta.index.str[0] == 'C']
meta.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all.csv")

# 算所有样本的所有亚类属于其大类的细胞比例
num2 = meta.groupby(['cluster', 'orig.ident']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.sub_cluster):
    x = meta[meta['sub_cluster'] == i]
    num1 = x.groupby('orig.ident').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)
t.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/percent_new/percent.csv")
```

# 按照差异基因选
```python
all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all.csv", index_col = 0)
meta = all[all.cluster == 'Mye']
mye_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/mye_res.csv", index_col=0, encoding = 'gbk')
marker = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/markers_batched2.csv")

sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.sub_cluster))
sg = sg.fillna(0)

sg = sg.replace(' ', '_')

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, mye_res.loc[mye_res.cluster == t.cluster].index] = score
sg = sg.loc[(sg.index.str[0:3] != 'RPL') & (sg.index.str[0:3] != 'RPS') & (sg.index.str[0:3] != 'MT-') & (sg.index.str[0:3] != 'RLS')]
for i in sg.columns:
    t = sg.loc[:, i]
    res = pd.DataFrame(t.loc[t > 0.3].sort_values(ascending = False)[0:30])
    res.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/crc_mye_score/"+ i + "_pct.csv")
```