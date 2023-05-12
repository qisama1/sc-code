tcga_meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/brca_meta.csv", index_col=0)
tpm = pd.read_csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/all-tpm_unstranded.csv", index_col=0).T
sampleid = pd.read_csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/sample_sheet.csv", index_col=0)
sampleid.index = sampleid['Sample ID']
# 去重
brca_counts = tpm.loc[sampleid.loc[sampleid['Case ID'].isin(tcga_meta.index)].index]

brca_counts = brca_counts[~brca_counts.index.duplicated()]

# 计算
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/all.csv", index_col=0)
meta.sub_cluster = meta.sub_cluster.str.replace(' ', '_').str.replace('-', '_')
s_celltype = pd.DataFrame(index = brca_counts.index, columns = set(meta.sub_cluster)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.sub_cluster == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/bc_mye_score/"+ i + "_pct.csv", index_col=0)
        cur = brca_counts.loc[:, brca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/TNK/bc_tnk_score/"+ i + "_pct.csv", index_col=0)
        cur = brca_counts.loc[:, brca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Fib/bc_fib_score/"+ i + "_pct.csv", index_col = 0)
        cur = brca_counts.loc[:, brca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Epi/bc_epi_score/"+ i + "_pct.csv", index_col = 0)
        cur = brca_counts.loc[:, brca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/bc_end_score/"+ i + "_pct.csv", index_col=0)
        cur = brca_counts.loc[:, brca_counts.columns.isin(res.index)]
        s_celltype.loc[:, i] = cur.mean(axis=1)

# 结合一下percent
s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/module/tpm_compare/s_celltype_gene.csv")

module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc_modules.csv", encoding = 'gbk')
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/bc3_t_cellpercent_sample.csv", index_col=0)
percent = percent.T
percent.columns = percent.columns.str.replace(' ', '_').str.replace('-', '_')

s_celltype.columns = s_celltype.columns.str.replace(' ', "_")

module = module.replace('Venous-1', 'Venous_1')

percent2 = percent.mean()
# s_celltype = s_celltype.replace('Naïve_B', 'Naive_B')
# percent2.index = percent2.index.str.replace('Naïve_B', 'Naive_B')

s_celltype.columns = s_celltype.columns.str.replace('anti_CAFs', 'ap_CAFs')
s_celltype.columns = s_celltype.columns.str.replace('Stem', 'Basal')


res = pd.DataFrame(index = s_celltype.index, columns = module.columns)
for p in res.index:
    for i in module.columns:
        s_subtme = 0
        for j in module.loc[:, i].dropna():
            s_subtme += (percent2[j] * s_celltype.loc[p, j]) ** 0.5
        res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())

tcga_meta['clinical_stage', 'neoplasm_histologic_grade']
# bcr_patient_barcode,gender,vital_status,                            
#                 days_to_death,days_to_last_followup,
#                 person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
#                 laterality,neoplasm_histologic_grade,pathologic_stage

## 区分癌症和正常
sampleid = sampleid[~sampleid.index.duplicated()]

res['type'] = sampleid.loc[res.index, 'Sample Type']
res.loc[res.type != 'Solid Tissue Normal', 'type'] = 'Tumor'
res.loc[res.type == 'Solid Tissue Normal', 'type'] = 'Normal'

res.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/module/tpm_compare/tpm_module_type.csv")


res = res.loc[res.type == 'Tumor']
meta_use = tcga_meta.loc[sampleid.loc[res.index, 'Case ID'], ['neoplasm_histologic_grade', 'pathologic_stage', 'person_neoplasm_cancer_status']]

meta_use = meta_use.fillna("nomsg")
meta_use = meta_use.loc[(meta_use.pathologic_stage != 'Stage X') & (meta_use.pathologic_stage != 'nomsg')]
meta_use['stage'] = 'Early'
meta_use.loc[meta_use.pathologic_stage.str.contains('Stage III') | meta_use.pathologic_stage.str.contains('Stage IV'), 'stage'] = 'Late'
meta_use['stage2'] = 'Early'
meta_use.loc[meta_use.pathologic_stage.str.contains('Stage II') | meta_use.pathologic_stage.str.contains('Stage IV'), 'stage2'] = 'Late'

t = sampleid.loc[res.index]
meta_use = meta_use[~meta_use.index.duplicated()]
res2 = res.loc[t[t['Case ID'].isin(meta_use.index)].index]
res2['stage'] = meta_use.loc[sampleid.loc[res2.index,'Case ID'], 'stage'].values
res2['stage2'] = meta_use.loc[sampleid.loc[res2.index,'Case ID'], 'stage2'].values

res2.to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/module/tpm_compare/tpm_module_stage.csv")