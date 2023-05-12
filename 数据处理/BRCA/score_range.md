# 取基因表达
```R
scRNA = qread("...")
write.csv(as.matrix(scRNA[gene$Gene]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/gene_df.csv")
```

# 算分数-用sg
```python
import pandas as pd
import os
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_scores/gene_df_filted.csv", index_col=0).T
meta.sub_cluster = meta.sub_cluster.str.replace(' ', '_').str.replace('-','_')
s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.sub_cluster)).fillna(0)

for i in range(30) :
    path = "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module_score_all/top" + i
    os.mkdir(path)
    for i in s_celltype.columns:
        major = meta.loc[meta.sub_cluster == i, 'cluster'].values[0]
        if (major == 'Mye'):
            res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/bc_mye_score/"+ i + "_pct.csv", index_col=0)
            res = res.iloc[0:i,]
            s_celltype.loc['score', i] = res.mean().values[0]
        if (major == 'TNK'):
            res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/TNK/bc_tnk_score/"+ i + "_pct.csv", index_col=0)
            res = res.iloc[0:i,]
            s_celltype.loc['score', i] = res.mean().values[0]
        if (major == 'Fib'):
            res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Fib/bc_fib_score/"+ i + "_pct.csv", index_col = 0)
            res = res.iloc[0:i,]
            s_celltype.loc['score', i] = res.mean().values[0]
        if (major == 'Epi'):
            res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Epi/bc_epi_score/"+ i + "_pct.csv", index_col = 0)
            res = res.iloc[0:i,]
            s_celltype.loc['score', i] = res.mean().values[0]
        if (major == 'End'):
            res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/bc_end_score/"+ i + "_pct.csv", index_col=0)
            res = res.iloc[0:i,]
            s_celltype.loc['score', i] = res.mean().values[0]
    s_celltype.to_csv(path + "s_celltype.csv")

    module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc_modules.csv", encoding = 'gbk')
    percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/percent_all.csv", index_col=0)
    module = module.replace('Venous-1', 'Venous_1')
    percent.columns = percent.columns.str.replace(' ', '_').str.replace('-', '_')
    percent.columns = percent.columns.str.replace('Stem', 'Basal')
    percent.columns = percent.columns.str.replace('anti_CAFs', 'ap_CAFs')

    s_celltype.columns = s_celltype.columns.str.replace(' ', "_")
    s_celltype.columns = s_celltype.columns.str.replace('anti_CAFs', 'ap_CAFs')
    s_celltype.columns = s_celltype.columns.str.replace('Stem', 'Basal')
    # 全部矫正一下
    module = module.replace(' ', '_')


    data = pd.DataFrame(index = percent.index, columns = percent.columns).fillna(0)
    for i in percent.index:
        for j in percent.columns:
            mod_col = j.replace(' ', '_')
            data.loc[i, mod_col] = percent.loc[i, j] * s_celltype.loc['score', j]

    res = pd.DataFrame(index = data.index, columns = module.columns)
    for p in res.index:
        for i in module.columns:
            s_subtme = 0
            for j in module.loc[:, i].dropna():
                s_subtme += (percent.loc[p, j].mean() * s_celltype.loc['score', j]) ** 0.5
            res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())

    res.to_csv(path + "s_res.csv")



p = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/dataset3_sampleinfo.csv", index_col=0)
```

# 用基因表达算
```python
s_celltype = pd.DataFrame(index = set(meta['orig.ident']), columns = set(meta.sub_cluster)).fillna(0)
gene_df = gene_df.loc[meta.index,]

for i in s_celltype.columns:
    major = meta.loc[meta.sub_cluster == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/bc_mye_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/TNK/bc_tnk_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Fib/bc_fib_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Epi/bc_epi_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/bc_end_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)

s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_gene_scores/s_celltype_gene.csv")

#module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/module.csv", encoding = 'gbk')
#percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/percent.csv", index_col=0)

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

res.loc[:, ['type', 'grade', 'stage']] = compare2.loc[:, ['type', 'grade', 'stage']]
```