tcga = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/tcga_tpm_log.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/pdac_sample_meta2.csv", index_col = 0).T
tcga_use = tcga.loc[meta.index,]
res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/tpm_module_type.csv", index_col=0)
an_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/4-12-tcga_heatmap-2/PDAC/an_gene.csv")
gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/ssgsea/pdac_ssgsva.csv", index_col = 0).T


gsva.index = gsva.index.str.replace('.', '-')
pathways = ['EMT', 'Immune escape',  'Pro-inflammatory', 'Matrix remodeling', 'Angiogenesis', 'Tumor proliferation', 'Hypoxia', 'Cellular?senescence']
corr_res = res.loc[:, ['subTME3']]
corr_res2 = gsva.loc[:, pathways]
corr_res2.columns = ['EMT', 'Immune_escape',  'Pro-inflammatory', 'Matrix remodeling', 'Angiogenesis', 'Tumor proliferation', 'Hypoxia', 'Cellular?senescence']
corr_res = pd.concat([corr_res, corr_res2], axis=1)

for pair in an_gene.gene:
    p = pd.DataFrame(index = corr_res.index, columns = [pair])
    if (len(pair.split('_')) <= 1) :
        gene1 = pair.split('|')[0]
        gene2 = pair.split('|')[1]
        p.loc[:, pair] = (tcga_use.loc[:, gene1] + tcga_use.loc[:, gene2]) / 2
    else:
        gene1 = pair.split('|')[0]
        gene2 = pair.split('|')[1].split('_')[0]
        gene3 = pair.split('|')[1].split('_')[1]
        p.loc[:, pair] = (tcga_use.loc[:, gene1] + tcga_use.loc[:, gene2] + tcga_use.loc[:, gene3]) / 3
    corr_res = pd.concat([corr_res, p], axis=1)
corr_res = corr_res.loc[corr_res.sort_values('subTME3',ascending=False).index,].T

import sklearn.preprocessing as preprocessing
basic_res = pd.DataFrame(preprocessing.scale(corr_res, axis=1), index = corr_res.index, columns=corr_res.columns)
basic_res[~basic_res.index.duplicated()].to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/4-12-tcga_heatmap-2/PDAC/heatmap_data.csv")
res.loc[:, ['type']].to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/4-12-tcga_heatmap-2/PDAC/an.csv")