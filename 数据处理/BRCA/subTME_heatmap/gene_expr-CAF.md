# gene expr

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/bc/an_pair.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/bc/gene_df.csv")
```

```py
# 1. 获取每行，然后判断类型，做一个简单的if/else即可，对于不同的类型对应不同的细胞类型的group
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/bc/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['sub_cluster']
data['sample'] = meta['orig.ident']


gene_used = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/bc/an_pair.csv")
basic_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_scores/module_compare/sg_diff_type.csv", index_col = 0)
basic_res = basic_res.loc[:, ['subTME1']]

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/t.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/mye.csv", index_col = 0).T
gsva_end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/end.csv", index_col = 0).T
gsva_fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/caf.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/epi.csv", index_col = 0).T


gsva_fib['sample'] = meta['orig.ident']
major = meta.loc[meta.cluster == 'Fib']
p = pd.DataFrame(gsva_fib.loc[gsva_fib.index.isin(major.index)].groupby('sample').mean().loc[:, ['myCAF']])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_all['sample'] = meta['orig.ident']
major = meta.loc[meta.cluster == 'Epi']
p = pd.DataFrame(gsva_all.loc[gsva_all.index.isin(major.index)].groupby('sample').mean().loc[:, ])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_epi['sample'] = meta['orig.ident']
major = meta.loc[meta.cluster == 'Epi']
p = pd.DataFrame(gsva_epi.loc[gsva_epi.index.isin(major.index)].groupby('sample').mean().loc[:, ['pEMT']])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_mye['sample'] = meta['orig.ident']
major = meta.loc[meta.sub_cluster.str.contains('Macro')]
p = pd.DataFrame(gsva_mye.loc[gsva_mye.index.isin(major.index)].groupby('sample').mean().loc[:, ])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_end['sample'] = meta['orig.ident']
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
basic_res = basic_res.loc[basic_res.sort_values('subTME1',ascending=False).index,].T
basic_res = basic_res[~basic_res.index.duplicated()]

idxs = ['subTME-MRM', 'myCAF', 'Proliferation', 'pEMT', 'angiogenic EC',
       'CAF_COL1A1', 'CAF_COL1A2', 'CAF_COL6A1', 'CAF_COL6A2', 'CAF_COL6A3',
       'CAF_FN1', 'CAF_THBS1', 'Epi_CD55', 'Epi_CXCL2', 'Epi_CXCL8',
       'Epi_CD47', 'Epi_CD44', 'Epi_IGF1R', 'Epi_CD24', 'Epi_NOTCH2',
       'Epi_SORL1', 'Epi_SLC1A5', 'Epi_TNFSF10', 'Epi_SDC1', 'Epi_SDC4',
       'TAM_NR3C1', 'TAM_SIRPA', 'TAM_TYROBP', 'TAM_ADGRE5', 'TAM_DPP4']
basic_res.index = idxs
import sklearn.preprocessing as preprocessing
res = pd.DataFrame(preprocessing.scale(basic_res, axis=1), index = basic_res.index, columns=basic_res.columns)
```