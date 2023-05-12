# 取基因表达
```R
scRNA = qread("...")
write.csv(as.matrix(scRNA@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/gene_df.csv")
```

# 算分数-用sg
```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/all_t_n.csv", index_col=0)
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/gene_df.csv", index_col=0).T
s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.ident)).fillna(0)
res_type = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_diff/s_celltype/res_type.csv", index_col = 0)
res_stage = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_diff/s_celltype/res_stage.csv", index_col = 0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/mye/19_mye_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/tnk/19_tnk_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/bcell/19_b_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/fib/19_fib_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/epi/19_epi_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/end/19_end_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/19/module/s_celltype.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/module.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/percent/percent.csv", index_col=0)

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
res.columns = res.columns.str.replace('module', 'subTME')
```

# 用基因表达算
```java
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/all_t_n.csv", index_col=0)
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/gene_df.csv", index_col=0).T
s_celltype = pd.DataFrame(index = set(meta['Sample']), columns = set(meta.ident)).fillna(0)
gene_df = gene_df.loc[meta.index,]
res_type = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_diff/s_celltype/res_type.csv", index_col = 0)
res_stage = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_diff/s_celltype/res_stage.csv", index_col = 0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/mye/19_mye_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'Sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/tnk/19_tnk_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'Sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/bcell/19_b_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'Sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/fib/19_fib_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'Sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/epi/19_epi_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'Sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/end/19_end_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.ident == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'Sample']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)

s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_gene_diff_range/s_celltype_gene.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/module/module.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/percent/percent.csv", index_col=0)

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
res.columns = res.columns.str.replace('module', 'subTME')
res.loc[res_type.index, 'type'] = res_type['type']
res.to_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_diff_range/" + str(step) + "_res_type.csv")
res = res.loc[res_stage.index, ]
res.loc[res_stage.index, ['stage1', 'stage2']] = res_stage.loc[:, ['stage1', 'stage2']]
res.to_csv("/public/home/yuwenqi/sc-data/selected/19/module/subTME_diff_range/" + str(step) + "_res_grade.csv")
```