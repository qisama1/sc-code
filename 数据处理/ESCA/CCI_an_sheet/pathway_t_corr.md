```python
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/t.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/epi.csv", index_col = 0).T
gsva_caf = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/caf.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/mye.csv", index_col = 0).T
gsva_end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/end.csv", index_col = 0).T

gene_corr = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gene_corr.csv")
# T
res = pd.DataFrame(index = gene_corr['T'].dropna(), columns = gsva_t.columns).fillna(0)
for gene in gene_corr['T'].dropna():
    for pathway in gsva_t.columns:
        meta_idx = meta_t[meta_t.ident.str.contains('CD8')].index
        idx = gsva_t.index[gsva_t.index.isin(meta_idx)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_t.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
idx = gsva_t.index[gsva_t.index.isin(meta_idx)]
gene_df.loc[idx, gene_corr['T'].dropna()].corr('spearman')
# Epi
res = pd.DataFrame(index = gene_corr['Epi'].dropna(), columns = gsva_epi.columns).fillna(0)
for gene in gene_corr['Epi'].dropna():
    for pathway in gsva_epi.columns:
        idx = gsva_epi.index[gsva_epi.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_epi.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/epi_corr.csv")
idx = gsva_epi.index[gsva_epi.index.isin(meta_t.index)]
gene_df.loc[idx, gene_corr['Epi'].dropna()].corr('spearman').to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/epi_gene_corr.csv")
# Mye
res = pd.DataFrame(index = gene_corr['Mye'].dropna(), columns = gsva_mye.columns).fillna(0)
for gene in gene_corr['Mye'].dropna():
    for pathway in gsva_mye.columns:
        idx = gsva_mye.index[gsva_mye.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_mye.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/mye_corr.csv")
idx = gsva_mye.index[gsva_mye.index.isin(meta_t.index)]
gene_df.loc[idx, gene_corr['Mye'].dropna()].corr('spearman').to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/mye_gene_corr.csv")
# CAF
res = pd.DataFrame(index = gene_corr['CAF'].dropna(), columns = gsva_caf.columns).fillna(0)
for gene in gene_corr['CAF'].dropna():
    for pathway in gsva_caf.columns:
        idx = gsva_caf.index[gsva_caf.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_caf.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/caf_corr.csv")
idx = gsva_caf.index[gsva_caf.index.isin(meta_t.index)]
gene_df.loc[idx, gene_corr['CAF'].dropna()].corr('spearman').to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/caf_gene_corr.csv")
# End
res = pd.DataFrame(index = gene_corr['End'].dropna(), columns = gsva_end.columns).fillna(0)
for gene in gene_corr['End'].dropna():
    for pathway in gsva_end.columns:
        idx = gsva_end.index[gsva_end.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_end.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/end_corr.csv")
idx = gsva_end.index[gsva_end.index.isin(meta_t.index)]
gene_df.loc[idx, gene_corr['End'].dropna()].corr('spearman').to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/end_gene_corr.csv")
# Epi_all
res = pd.DataFrame(index = gene_corr['Epi'].dropna(), columns = gsva_all.columns).fillna(0)
for gene in gene_corr['Epi'].dropna():
    for pathway in gsva_all.columns:
        idx = gsva_all.index[gsva_all.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_all.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/epi_proliferation_corr.csv")
idx = gsva_all.index[gsva_all.index.isin(meta_t.index)]
gene_df.loc[idx, gene_corr['Epi'].dropna()].corr('spearman').to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/epi_proliferation_gene_corr.csv")
```