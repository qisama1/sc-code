# gene expr

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/esca/an_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/esca/gene_df.csv")
```

```py
# 1. 获取每行，然后判断类型，做一个简单的if/else即可，对于不同的类型对应不同的细胞类型的group
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/esca/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['ident']
data['sample'] = meta['sample']

gene_used = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/esca/an_gene.csv")
basic_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_compare.csv", index_col = 0)
basic_res = basic_res.loc[:, ['subTME2']]

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/t.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/mye.csv", index_col = 0).T
gsva_fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/caf.csv", index_col = 0).T
gsva_end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/end.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/epi.csv", index_col = 0).T

gsva_fib['sample'] = meta['sample']
major = meta.loc[meta.cluster == 'Fib']
p = pd.DataFrame(gsva_fib.loc[gsva_fib.index.isin(major.index)].groupby('sample').mean().loc[:, ['myCAF']])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_all['sample'] = meta['sample']
major = meta.loc[meta.cluster == 'Epi']
p = pd.DataFrame(gsva_all.loc[gsva_all.index.isin(major.index)].groupby('sample').mean().loc[:, ])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_epi['sample'] = meta['sample']
major = meta.loc[meta.cluster == 'Epi']
p = pd.DataFrame(gsva_epi.loc[gsva_epi.index.isin(major.index)].groupby('sample').mean().loc[:, ['pEMT']])
basic_res = pd.concat([basic_res, p], axis=1)

# gsva_mye['sample'] = meta['sample']
# major = meta.loc[meta.sub_cluster.str.contains('Macro')]
# p = pd.DataFrame(gsva_mye.loc[gsva_mye.index.isin(major.index)].groupby('sample').mean().loc[:, ])
# basic_res = pd.concat([basic_res, p], axis=1)

gsva_end['sample'] = meta['sample']
major = meta.loc[meta.cluster == 'End']
p = pd.DataFrame(gsva_end.loc[gsva_end.index.isin(major.index)].groupby('sample').mean().loc[:, ['angiogenic EC']])
basic_res = pd.concat([basic_res, p], axis=1)

for i in gene_used.index:
    v = gene_used.iloc[i, ]
    g = v.gene
    if (g not in data.columns):
        continue
    if (v.type == 'TAM'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('Macro')].index,]
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        p.columns = ["TAM_" + g]
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'CD8+T'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('CD8')].index,]
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        p.columns = ["CD8_" + g]
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'Treg'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('Treg')].index,]
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        p.columns = ["Treg_" + g]
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'Epi'):
        major = meta.loc[meta.cluster == 'Epi']
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        p.columns = ["Epi_" + g]
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'CAF'):
        major = meta.loc[meta.cluster == 'Fib']
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        p.columns = ["CAF_" + g]
        basic_res = pd.concat([basic_res, p], axis=1)
basic_res = basic_res.fillna(0)
basic_res = basic_res.loc[basic_res.sort_values('subTME2',ascending=False).index,].T

idxs =  ['subTME-MRM', 'myCAF', 'Proliferation', 'pEMT', 'angiogenic EC',
       'CAF_COL1A1', 'CAF_COL1A2', 'CAF_COL6A1', 'CAF_COL6A2', 'CAF_COL6A3',
       'CAF_FN1', 'CAF_LAMA4', 'CAF_LAMB1', 'CAF_LAMC1', 'CAF_NAMPT',
       'CAF_THBS1', 'CAF_THBS2', 'CAF_THY1', 'Epi_ITGA2', 'Epi_ITGA3',
       'Epi_ITGB1', 'Epi_SDC1', 'Epi_SDC4']

basic_res.index = idxs
import sklearn.preprocessing as preprocessing
res = pd.DataFrame(preprocessing.scale(basic_res, axis=1), index = basic_res.index, columns=basic_res.columns)
pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/subtem_compare.csv", index_col = 0).loc[:, ['type']].to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/esca/an.csv")
```