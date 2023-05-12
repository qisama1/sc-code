# 取基因表达
```R
scRNA = qread("...")
write.csv(as.matrix(scRNA@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/gene_df.csv")
```

# 算分数-用sg
```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all2.csv", index_col=0)
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/gene_df.csv", index_col=0).T
s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.ident)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/mye/9_mye_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/tnk/9_tnk_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/bcell/9_b_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/fib/9_fib_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/epi/9_epi_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/endo/9_end_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/s_celltype.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/module.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/percent2.csv", index_col=0)

data = pd.DataFrame(index = percent.index, columns = percent.columns).fillna(0)
for i in data.index:
    for j in data.columns:
        data.loc[i, j] = percent.loc[i, j] * s_celltype.loc['score', j]

res = pd.DataFrame(index = data.index, columns = module.columns)
for p in res.index:
    for i in module.columns:
        s_subtme = 0
        for j in module.loc[:, i].dropna():
            s_subtme += (percent.loc[p, j].mean() * s_celltype.loc['score', j]) ** 0.5
        res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())
```

# 用基因表达算
```java
s_celltype = pd.DataFrame(index = set(meta['sample']), columns = set(meta.ident)).fillna(0)
gene_df = gene_df.loc[meta.index,]

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/mye/9_mye_score/"+ i + "_pct.csv", index_col=0)
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/tnk/9_tnk_score/"+ i + "_pct.csv", index_col=0)
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/bcell/9_b_score/"+ i + "_pct.csv", index_col = 0)
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/fib/9_fib_score/"+ i + "_pct.csv", index_col = 0)
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/epi/9_epi_score/"+ i + "_pct.csv", index_col = 0)
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/endo/9_end_score/"+ i + "_pct.csv", index_col=0)
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)

s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/s_celltype_gene.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/module.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/percent2.csv", index_col=0)

data = pd.DataFrame(index = percent.index, columns = percent.columns).fillna(0)
for i in data.index:
    for j in data.columns:
        data.loc[i, j] = percent.loc[i, j] * s_celltype.loc[i, j]

res = pd.DataFrame(index = data.index, columns = module.columns)
for p in res.index:
    for i in module.columns:
        s_subtme = 0
        for j in module.loc[:, i].dropna():
            s_subtme += (percent.loc[p, j].mean() * s_celltype.loc[p, j]) ** 0.5
        res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())
compare2 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_compare2.csv", index_col=0)
res.loc[:, ['type', 'grade', 'stage']] = compare2.loc[:, ['type', 'grade', 'stage']]
```