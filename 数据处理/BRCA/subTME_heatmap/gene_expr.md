# gene expr

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/subTME_heatmap/BC/an_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/subTME_heatmap/BC/gene_df.csv")


```

```py
# 1. 获取每行，然后判断类型，做一个简单的if/else即可，对于不同的类型对应不同的细胞类型的group
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/subTME_heatmap/BC/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['sub_cluster']
data['sample'] = meta['orig.ident']


gene_used = pd.read_csv("/public/home/yuwenqi/sc-data/selected/subTME_heatmap/BC/an_needed.csv")
basic_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_scores/module_compare/sg_diff_type.csv", index_col = 0)
basic_res = basic_res.loc[:, ['subTME2', 'subTME1']]

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/t.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/mye.csv", index_col = 0).T



gsva_mye['sample'] = meta['orig.ident']
major = meta.loc[meta.sub_cluster.str.contains('Macro')]
p = pd.DataFrame(gsva_mye.loc[gsva_mye.index.isin(major.index)].groupby('sample').mean().loc[:, ])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_t['sample'] = meta['orig.ident']
major = meta.loc[meta.sub_cluster.str.contains('CD8')]
p = pd.DataFrame(gsva_t.loc[gsva_t.index.isin(major.index)].groupby('sample').mean().loc[:, 'Exhausted'])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_all['sample'] = meta['orig.ident']
major = meta.loc[meta.cluster == 'Epi']
p = pd.DataFrame(gsva_all.loc[gsva_all.index.isin(major.index)].groupby('sample').mean().loc[:, ])
basic_res = pd.concat([basic_res, p], axis=1)

for i in gene_used.index:
    v = gene_used.iloc[i, ]
    g = v.gene
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
basic_res = basic_res[~basic_res.index.duplicated()]

idxs = ['subTME2', 'subTME1', 'M1_Macrophage_Polarization', 'M2_Macrophage_Polarization',
       'Pro_inflammatory_Mye', 'Anti-inflammatory_Mye', 'Exhausted', 'Proliferation', 'TAM_TNF',
       'TAM_CD47', 'TAM_P2RY6', 'TAM_PLAUR', 'TAM_IL10', 'TAM_CCL7',
       'TAM_CD86', 'TAM_SIGLEC10', 'TAM_SIGLEC1', 'TAM_CCL2', 'TAM_NECTIN2',
       'TAM_MRC1', 'TAM_CXCL9', 'TAM_CXCL8', 'TAM_GRN', 'TAM_NR3C1',
       'TAM_SIRPA', 'TAM_TYROBP', 'TAM_ADGRE5', 'TAM_DPP4', 'TAM_LRP1',
       'CD8_CD52', 'CD8_SPN', 'CD8_NR3C1', 'CD8_TIGIT', 'CD8_PTPRC',
       'CD8_CXCR3', 'CD8_TNFRSF1B', 'Treg_NR3C1', 'Treg_CTLA4', 'Treg_ICOS',
       'Treg_SIRPG', 'Epi_CD55', 'Epi_CXCL2', 'Epi_MDK', 'Epi_CXCL8',
       'Epi_CD47', 'Epi_CD44', 'Epi_EGFR', 'Epi_SDC1', 'Epi_ITGA2',
       'Epi_ITGA3', 'Epi_ITGB1', 'Epi_IGFBP3', 'Epi_RPS19', 'Epi_SLC7A1',
       'Epi_IDE', 'Epi_FGFR2', 'Epi_ANXA1', 'Epi_PLAU', 'CAF_NAMPT',
       'CAF_PLAU', 'CAF_THBS1', 'CAF_COL1A1', 'CAF_COL1A2', 'CAF_COL6A3',
       'CAF_COL6A2', 'CAF_THY1', 'CAF_LAMC1', 'CAF_LAMB1', 'CAF_IGFBP3']
basic_res.index = idxs
import sklearn.preprocessing as preprocessing
res = pd.DataFrame(preprocessing.scale(basic_res, axis=1), index = basic_res.index, columns=basic_res.columns)
```