# get expr
```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/needed_gene.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/gene_df.csv")
```

# corr
```py
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/gene_df.csv", index_col=0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)

meta_t = meta.loc[meta.Tissue == 'Tumor']

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/t.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/epi.csv", index_col = 0).T
gsva_caf = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/caf.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/mye.csv", index_col = 0).T


def spearmanr_pval(x,y):
    if (spearmanr(x, y)[1] > 0.05):
        return 0
    else :
        return spearmanr(x, y)[0]

# T
an_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/an_gene_t.csv")
res = pd.DataFrame(index = an_t.gene.dropna(), columns = gsva_t.columns).fillna(0)
meta_idx = meta_t[meta_t.sub_cluster.str.contains('CD8')].index
for gene in an_t.gene.dropna():
    for pathway in gsva_t.columns:
        meta_idx = meta_t[meta_t.sub_cluster.str.contains('CD8')].index
        idx = gsva_t.index[gsva_t.index.isin(meta_idx)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_t.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]

idx = gsva_t.index[gsva_t.index.isin(meta_idx)]
gene_df.loc[idx, an_t.gene.dropna()].corr(method = spearmanr_pval)
# Epi
an_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/an_gene_epi.csv")
res = pd.DataFrame(index = an_epi.gene.dropna(), columns = gsva_epi.columns).fillna(0)
for gene in an_epi.gene.dropna():
    if (gene not in gene_df.columns):
        continue
    for pathway in gsva_epi.columns:
        idx = gsva_epi.index[gsva_epi.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_epi.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/corr_epi.csv")
idx = gsva_epi.index[gsva_epi.index.isin(meta_t.index)]
gene_df.loc[idx, an_t.gene.dropna()[an_t.gene.dropna().isin(gene_df.columns)]].corr(method = spearmanr_pval).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/corr_epi_gene.csv")
# Mye
an_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/an_gene_tam.csv")
res = pd.DataFrame(index = an_mye.gene.dropna(), columns = gsva_mye.columns).fillna(0)
meta_idx = meta_t[meta_t.sub_cluster.str.contains('Macro')].index
for gene in an_mye.gene.dropna():
    if (gene not in gene_df.columns):
        continue
    for pathway in gsva_mye.columns:
        idx = gsva_mye.index[gsva_mye.index.isin(meta_idx)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_mye.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/corr_mye.csv")
idx = gsva_mye.index[gsva_mye.index.isin(meta_idx)]
gene_df.loc[idx, an_mye.gene.dropna()[an_mye.gene.dropna().isin(gene_df.columns)]].corr(method = spearmanr_pval).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/corr_mye_gene.csv")

# Epi_all
res = pd.DataFrame(index = an_epi.gene.dropna(), columns = gsva_all.columns).fillna(0)
for gene in an_epi.gene.dropna():
    if (gene not in gene_df.columns):
        continue
    for pathway in gsva_all.columns:
        idx = gsva_all.index[gsva_all.index.isin(meta_t.index)]
        corr = spearmanr(gene_df.loc[idx, gene], gsva_all.loc[idx, pathway])
        if (corr[1] < 0.05):
            res.loc[gene, pathway] = corr[0]
            
res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/corr_epi_proliferation.csv")
idx = gsva_all.index[gsva_all.index.isin(meta_t.index)]
gene_df.loc[idx, an_t.gene.dropna()[an_t.gene.dropna().isin(gene_df.columns)]].corr(method = spearmanr_pval).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/BC/corr_epi_proliferation_gene.csv")

# gene_corr
an_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/an_gene_gene.csv")
res = pd.DataFrame(index = an_gene.gene_1.dropna(), columns = an_gene.gene_2.dropna()).fillna(0)
meta_idx = meta_t[meta_t.sub_cluster.str.contains('Macro')].index
for i in res.index:
    if (i not in gene_df.columns):
        continue
    for j in res.columns:
        if (j not in gene_df.columns):
            continue
        r = spearmanr(gene_df.loc[meta_idx, i], gene_df.loc[meta_idx, j])
        if (r[1] < 0.05):
            res.loc[i, j] = r[0]
```