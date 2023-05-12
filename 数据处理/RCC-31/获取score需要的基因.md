# py
```py
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/all_t_n.csv", index_col=0)
genes = set()
s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.ident)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/mye/31_mye_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/tnk/31_tnk_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/bcell/31_b_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/fib/31_fib_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/epi/31_epi_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/31/end/31_end_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
pd.DataFrame(genes, columns = ['gene']).to_csv("/public/home/yuwenqi/sc-data/selected/31/module/gene_neeeded.csv")

```

# R
```R
gene = read.csv("/public/home/yuwenqi/sc-data/selected/31/module/gene_neeeded.csv")
scRNA = qread("/public/home/yuwenqi/sc-data/selected/31/all.qs")
write.csv(as.matrix(scRNA[gene$gene]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/31/module/gene_df.csv")
```