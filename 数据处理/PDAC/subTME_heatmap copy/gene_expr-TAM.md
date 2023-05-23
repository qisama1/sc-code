# gene expr

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME-heatmap-TAM-cancer/pdac/an_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/subTME-heatmap-TAM-cancer/pdac/gene_df.csv")
```

```py
# 1. 获取每行，然后判断类型，做一个简单的if/else即可，对于不同的类型对应不同的细胞类型的group
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME-heatmap-TAM-cancer/pdac/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
data = data.loc[meta.index, ]
data['sub_cluster'] = meta['ident']
data['sample'] = meta['sample']

gene_used = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME-heatmap-TAM-cancer/pdac/an_gene.csv")
basic_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cluster_sg_scores/s_celltype_gene_type.csv", index_col = 0)
basic_res = basic_res.loc[:, ['subTME3']]

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/t.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/mye.csv", index_col = 0).T
gsva_fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/caf.csv", index_col = 0).T
gsva_end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/end.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/epi.csv", index_col = 0).T

gsva_mye['sample'] = meta['sample']
major = meta.loc[meta.ident.str.contains('Macro')]
p = pd.DataFrame(gsva_mye.loc[gsva_mye.index.isin(major.index)].groupby('sample').mean().loc[:, ['M2 Macrophage Polarization']])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_mye['sample'] = meta['sample']
major = meta.loc[meta.ident.str.contains('Macro')]
p = pd.DataFrame(gsva_mye.loc[gsva_mye.index.isin(major.index)].groupby('sample').mean().loc[:, ['Anti-inflammatory (Mye)']])
basic_res = pd.concat([basic_res, p], axis=1)

gsva_t['sample'] = meta['sample']
major = meta.loc[meta.ident.str.contains('CD8')]
p = pd.DataFrame(gsva_t.loc[gsva_t.index.isin(major.index)].groupby('sample').mean().loc[:, ['Exhausted']])
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
basic_res = basic_res.loc[basic_res.sort_values('subTME3',ascending=False).index,].T

idxs =  ['subTME-IS', 'M2 Macrophage Polarization', 'Anti-inflammatory (Mye)',
       'Exhausted', 'Epi_ANXA1', 'Epi_IDE', 'Epi_MDK', 'Epi_PLAU', 'Epi_RPS19',
       'Epi_SLC7A1', 'TAM_C3AR1', 'TAM_C5AR1', 'TAM_CCL3', 'TAM_CCL4',
       'TAM_CD40', 'TAM_CD47', 'TAM_CXCL8', 'TAM_FPR3', 'TAM_LRP1', 'TAM_MRC1',
       'TAM_NECTIN2', 'TAM_PLAUR', 'TAM_PTPRC', 'TAM_SPP1', 'TAM_SIRPA',
       'CD8_CSF1', 'CD8_NR3C1', 'CD8_SIRPA', 'CD8_TIGIT']

basic_res.index = idxs
import sklearn.preprocessing as preprocessing
res = pd.DataFrame(preprocessing.scale(basic_res, axis=1), index = basic_res.index, columns=basic_res.columns)
pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cluster_sg_scores/s_celltype_gene_type.csv", index_col = 0).loc[:, ['type']].to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/subTME-heatmap-TAM-cancer/pdac/an.csv")
```