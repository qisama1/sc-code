import pandas as pd

gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/pdac/gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
genes = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/cluster_genes_all_used.csv")

gene_df = gene_df.loc[meta.index, genes.gene[genes.gene.isin(gene_df.columns)]]
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/module.csv")
res = pd.DataFrame()

for m in ['module1', 'module2', 'module3', 'module4']:
    for cluster in module[m].dropna():
        res = pd.concat([res, pd.DataFrame([[cluster, m]], index = [cluster], columns = ['sub_cluster', 'subTME'])])

for pathway in gene_df.columns:
    for idx in res.index:
        cur = meta.loc[meta.ident == idx]
        res.loc[idx, pathway] = gene_df.loc[cur.index, pathway].mean()
res = res.iloc[:, 1:]

res = res.replace('module1', 'subTME-PSE')
res = res.replace('module2', 'subTME-PSE')
res = res.replace('module3', 'subTME-IS+MRM')
res = res.replace('module4', 'subTME-MRM')

res.to_csv("/public/home/yuwenqi/sc-data/selected/review/pdac/subTME_diff_plot/plot_data.csv")

