# gene expr

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_used.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/subTME_heatmap/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/subTME_heatmap/CRC/gene_df.csv")
```

```py
# 1. 获取每行，然后判断类型，做一个简单的if/else即可，对于不同的类型对应不同的细胞类型的group
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/subTME_heatmap/CRC/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all3.csv", index_col=0)
meta.loc[meta.index.str.split('_').str[1] == 'N', 'orig.ident'] = meta['orig.ident'] + '_n'
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['sub_cluster']
data['sample'] = meta['orig.ident']


res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/cluster_sg_gene_score/s_celltype_res_type.csv", index_col = 0)
grade = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/cluster_sg_gene_score/s_celltype_res_grade.csv", index_col=0)

gene_used = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME_heatmap/crc/an_pair.csv")
basic_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/cluster_sg_scores/s_res_type.csv", index_col = 0)
basic_res = basic_res.loc[:, ['subTME1']]

gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/t.csv", index_col = 0).T
gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/all.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/mye.csv", index_col = 0).T
gsva_fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/caf.csv", index_col = 0).T

gsva_fib['sample'] = meta['orig.ident']
major = meta.loc[meta.cluster == 'Fib']
p = pd.DataFrame(gsva_fib.loc[gsva_fib.index.isin(major.index)].groupby('sample').mean().loc[:, ['myCAF']])
basic_res = pd.concat([basic_res, p], axis=1)

# gsva_mye['sample'] = meta['orig.ident']
# major = meta.loc[meta.sub_cluster.str.contains('Macro')]
# p = pd.DataFrame(gsva_mye.loc[gsva_mye.index.isin(major.index)].groupby('sample').mean().loc[:, ])
# basic_res = pd.concat([basic_res, p], axis=1)

gsva_all['sample'] = meta['orig.ident']
major = meta.loc[meta.cluster == 'Epi']
p = pd.DataFrame(gsva_all.loc[gsva_all.index.isin(major.index)].groupby('sample').mean().loc[:, ])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_t['sample'] = meta['orig.ident']
major = meta.loc[meta.sub_cluster.str.contains('CD8')]
p = pd.DataFrame(gsva_t.loc[gsva_t.index.isin(major.index)].groupby('sample').mean().loc[:, 'Exhausted'])
basic_res = pd.concat([basic_res, p], axis=1)



for i in gene_used.index:
    v = gene_used.iloc[i, ]
    g = v.gene
    if (v.type == 'TAM'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('Macro')].index,]
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'CD8+T'):
        major = meta.loc[data.loc[data.sub_cluster.str.contains('CD8')].index,]
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'Epi'):
        major = meta.loc[meta.cluster == 'Epi']
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        basic_res = pd.concat([basic_res, p], axis=1)
    if (v.type == 'CAF'):
        major = meta.loc[meta.cluster == 'Fib']
        p = pd.DataFrame(data.loc[data.index.isin(major.index)].groupby('sample').mean().loc[:, g])
        basic_res = pd.concat([basic_res, p], axis=1)
basic_res = basic_res.fillna(0)
basic_res = basic_res.loc[basic_res.sort_values('subTME_1',ascending=False).index,].T

idxs =  ['subTME_1', 'M1_Macrophage_Polarization', 'M2_Macrophage_Polarization',
       'Pro_inflammatory_Mye', 'Anti-inflammatory (Mye)', 'Exhausted',
       'Proliferation', 'ITGAM', 'ITGAX', 'ITGA6', 'INSR', 'CD36', 'TMEM219',
       'TNFSF13B', 'ARF1', 'MRC1', 'ITGB2', 'ITGB1', 'THBS1', 'CD40', 'PLD2',
       'PTPRC', 'LRP1', 'FPR1', 'FPR3', 'PLAUR', 'NECTIN2', 'CXCL8', 'SPP1',
       'C5AR1', 'CCL4', 'CCL3', 'CD47', 'TIMP1', 'C5AR1.1', 'TMEM219.1',
       'PTGER4', 'SIRPA', 'RPS19', 'CSF1', 'TIGIT', 'NR3C1', 'SDC1', 'ITGA2',
       'ITGA3', 'ITGB1.1', 'ITGB1.2', 'IGFBP3', 'RPS19.1', 'SLC7A1', 'IDE',
       'FGFR2', 'RPS19.2', 'MDK', 'ANXA1', 'PLAU', 'THBS1.1', 'COL1A1',
       'COL1A2', 'COL6A3', 'COL6A2', 'COL1A1.1', 'COL1A1.2', 'THY1', 'LAMC1',
       'LAMB1', 'NAMPT', 'IGFBP3.1']

basic_res.index = idxs
import sklearn.preprocessing as preprocessing
res = pd.DataFrame(preprocessing.scale(basic_res, axis=1), index = basic_res.index, columns=basic_res.columns)
```