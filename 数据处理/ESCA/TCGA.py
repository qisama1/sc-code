tcga_meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ESCA_meta.csv", index_col=0)
tpm = pd.read_csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/all-tpm_unstranded.csv", index_col=0).T
sampleid = pd.read_csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/sample_sheet.csv", index_col=0)

sampleid.index = sampleid['Sample ID']
esca_counts = tpm.loc[sampleid.loc[sampleid['Case ID'].isin(tcga_meta.index)].index]

## 区分正常和癌症
set(sampleid.loc[sampleid['Sample ID'].isin(esca_counts.index), 'Sample Type'])
# 去重
esca_counts = esca_counts[~esca_counts.index.duplicated()]

# 计算
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all2.csv", index_col=0)
s_celltype = pd.DataFrame(index = esca_counts.index, columns = set(meta.ident)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/mye/9_mye_score/"+ i + "_pct.csv", index_col=0)
        cur = esca_counts.loc[:, esca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/tnk/9_tnk_score/"+ i + "_pct.csv", index_col=0)
        cur = esca_counts.loc[:, esca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/bcell/9_b_score/"+ i + "_pct.csv", index_col = 0)
        cur = esca_counts.loc[:, esca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/fib/9_fib_score/"+ i + "_pct.csv", index_col = 0)
        cur = esca_counts.loc[:, esca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/epi/9_epi_score/"+ i + "_pct.csv", index_col = 0)
        cur = esca_counts.loc[:, esca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/endo/9_end_score/"+ i + "_pct.csv", index_col=0)
        cur = esca_counts.loc[:, esca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)

# 结合一下percent
s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/s_celltype_gene.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/module.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/percent2.csv", index_col=0)

percent2 = percent.mean()

res = pd.DataFrame(index = esca_counts.index, columns = module.columns)
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

res.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_type.csv")

## 区分分期
res = res.loc[res.type == 'Tumor']
meta_use = tcga_meta.loc[sampleid.loc[res.index, 'Case ID'], ['neoplasm_histologic_grade', 'pathologic_stage', 'person_neoplasm_cancer_status']]

meta_use = meta_use.fillna("nomsg")
meta_use = meta_use.loc[(meta_use.pathologic_stage != 'Stage X') & (meta_use.pathologic_stage != 'nomsg')]
meta_use['stage'] = 'Early'
meta_use.loc[meta_use.pathologic_stage.str.contains('Stage III') | meta_use.pathologic_stage.str.contains('Stage IV'), 'stage'] = 'Late'
meta_use['stage2'] = 'Early'
meta_use.loc[meta_use.pathologic_stage.str.contains('Stage II') | meta_use.pathologic_stage.str.contains('Stage IV'), 'stage2'] = 'Late'

res.index = sampleid.loc[res.index, 'Case ID']
res = res.loc[meta_use.index, ]
res.loc[:, ['stage', 'stage2']] = meta_use.loc[res.index, ['stage', 'stage2']]

res.drop('type', axis=1).to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_stage.csv")

meta_use = tcga_meta.loc[sampleid.loc[res.index, 'Case ID'], ['neoplasm_histologic_grade', 'pathologic_stage', 'person_neoplasm_cancer_status']]
meta_use = meta_use.loc[(meta_use.neoplasm_histologic_grade != 'GX') & (meta_use.neoplasm_histologic_grade != 'nomsg')]
meta_use['grade'] = 'Low'
meta_use.loc[meta_use.neoplasm_histologic_grade.str.contains('G3') | meta_use.neoplasm_histologic_grade.str.contains('G4'), 'grade'] = 'High'
meta_use['grade2'] = 'Low'
meta_use.loc[meta_use.neoplasm_histologic_grade.str.contains('G2') | meta_use.neoplasm_histologic_grade.str.contains('G3'), 'grade2'] = 'High'
res = res.loc[meta_use.index, ]
res.loc[:, ['grade', 'grade2']] = meta_use.loc[res.index, ['grade', 'grade2']]

res.drop('stage', axis=1).drop('stage2', axis=1).to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_grade.csv")