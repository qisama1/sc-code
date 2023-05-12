# py
```py
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/all_t_n.csv", index_col=0)
genes = set()
s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.ident)).fillna(0)

for i in s_celltype.columns:
    print(i)
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/mye/19_mye_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    elif (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/tnk/19_tnk_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    elif (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/bcell/19_b_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    elif (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/fib/19_fib_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    elif (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/epi/19_epi_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
    elif (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/end/19_end_score/"+ i + "_pct.csv", index_col=0)
        for gene in res.index:
            genes.add(gene)
pd.DataFrame(genes, columns = ['gene']).to_csv("/public/home/yuwenqi/sc-data/selected/19/module/gene_neeeded.csv")
```

# R
```R
gene = read.csv("/public/home/yuwenqi/sc-data/selected/19/module/gene_neeeded.csv")
scRNA = qread("/public/home/yuwenqi/sc-data/selected/19/all.qs")
write.csv(as.matrix(scRNA[gene$gene]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/19/module/gene_df.csv")
```