tcga_meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/pdac_meta.csv", index_col=0)
tpm = pd.read_csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/all-tpm_unstranded.csv", index_col=0).T
sampleid = pd.read_csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/sample_sheet.csv", index_col=0)
sampleid.index = sampleid['Sample ID']
pdac_counts = tpm.loc[sampleid.loc[sampleid['Case ID'].isin(tcga_meta.index)].index]

# 去重
pdac_counts = pdac_counts[~pdac_counts.index.duplicated()]

# 计算
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all.csv", index_col=0)
meta.ident = meta.ident.str.replace(' ', '_')
s_celltype = pd.DataFrame(index = pdac_counts.index, columns = set(meta.ident)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/mye/35_mye_score/"+ i + "_pct.csv", index_col=0)
        cur = pdac_counts.loc[:, pdac_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/tnk/35_tnk_score/"+ i + "_pct.csv", index_col=0)
        cur = pdac_counts.loc[:, pdac_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/bcell/35_b_score/"+ i + "_pct.csv", index_col = 0)
        cur = pdac_counts.loc[:, pdac_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/fib/35_fib_score/"+ i + "_pct.csv", index_col = 0)
        cur = pdac_counts.loc[:, pdac_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/epi/35_epi_score/"+ i + "_pct.csv", index_col = 0)
        cur = pdac_counts.loc[:, pdac_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/endo/35_end_score/"+ i + "_pct.csv", index_col=0)
        cur = pdac_counts.loc[:, pdac_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)

# 结合一下percent
s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/s_celltype_gene.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/module.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/percent.csv", index_col=0)
percent.columns = percent.columns.str.replace(' ', '_')
s_celltype.columns = s_celltype.columns.str.replace(' ', "_")
module = module.replace('Hypoxia?CAF', 'Hypoxia-CAF')
module = module.replace('Langerhans?like', 'Langerhans-like')

module = module.replace(' ', '_')
percent2 = percent.mean()
s_celltype = s_celltype.replace('Naïve_B', 'Naive_B')
percent2.index = percent2.index.str.replace('Naïve_B', 'Naive_B')

res = pd.DataFrame(index = s_celltype.index, columns = module.columns)
for p in res.index:
    for i in module.columns:
        s_subtme = 0
        for j in module.loc[:, i].dropna():
            s_subtme += (percent2[j] * s_celltype.loc[p, j]) ** 0.5
        res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())

## 区分癌症和正常
res['type'] = sampleid.loc[res.index, 'Sample Type']
res.loc[res.type != 'Solid Tissue Normal', 'type'] = 'Tumor'
res.loc[res.type == 'Solid Tissue Normal', 'type'] = 'Normal'

res.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/tpm_module_type.csv")

## 結果
res = res.loc[res.type == 'Tumor']
meta_use = tcga_meta.loc[sampleid.loc[res.index, 'Case ID'], ['neoplasm_histologic_grade', 'pathologic_stage', 'person_neoplasm_cancer_status']]
meta_use = meta_use[~meta_use.index.duplicated()]
meta_use = meta_use.fillna("nomsg")
meta_use = meta_use.loc[(meta_use.pathologic_stage != 'Stage X') & (meta_use.pathologic_stage != 'nomsg')]
meta_use['stage'] = 'Early'
meta_use.loc[meta_use.pathologic_stage.str.contains('Stage III') | meta_use.pathologic_stage.str.contains('Stage IV'), 'stage'] = 'Late'
meta_use['stage2'] = 'Early'
meta_use.loc[meta_use.pathologic_stage.str.contains('Stage II') | meta_use.pathologic_stage.str.contains('Stage IV'), 'stage2'] = 'Late'

t = sampleid.loc[res.index]
res2 = res.loc[t[t['Case ID'].isin(meta_use.index)].index]
res2['stage'] = meta_use.loc[sampleid.loc[res2.index,'Case ID'], 'stage'].values
res2['stage2'] = meta_use.loc[sampleid.loc[res2.index,'Case ID'], 'stage2'].values


res2.drop('type', axis=1).to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/tpm_module_stage.csv")

meta_use = tcga_meta.loc[sampleid.loc[res.index, 'Case ID'], ['neoplasm_histologic_grade', 'pathologic_stage', 'person_neoplasm_cancer_status']]
meta_use = meta_use.loc[(meta_use.neoplasm_histologic_grade != 'GX') & (meta_use.neoplasm_histologic_grade != 'nomsg')]
meta_use['grade'] = 'Low'
meta_use.loc[meta_use.neoplasm_histologic_grade.str.contains('G3') | meta_use.neoplasm_histologic_grade.str.contains('G4'), 'grade'] = 'High'
meta_use['grade2'] = 'Low'
meta_use.loc[meta_use.neoplasm_histologic_grade.str.contains('G2') | meta_use.neoplasm_histologic_grade.str.contains('G3'), 'grade2'] = 'High'

t = sampleid.loc[res.index]
meta_use = meta_use[~meta_use.index.duplicated()]
res2 = res.loc[t[t['Case ID'].isin(meta_use.index)].index]
res2['grade'] = meta_use.loc[sampleid.loc[res2.index,'Case ID'], 'grade'].values
res2['grade2'] = meta_use.loc[sampleid.loc[res2.index,'Case ID'], 'grade2'].values


res2.drop('type', axis=1).to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/tpm_module_grade.csv")