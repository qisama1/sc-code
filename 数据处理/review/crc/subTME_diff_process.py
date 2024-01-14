import pandas as pd

gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/crc/gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all3.csv", index_col = 0)
genes = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/cluster_genes_all_used.csv")

gene_df = gene_df.loc[meta.index, genes.gene[genes.gene.isin(gene_df.columns)]]
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module.csv", index_col=0)
res = pd.DataFrame()

for m in ['module1', 'module2', 'module3']:
    for cluster in module[m].dropna():
        res = pd.concat([res, pd.DataFrame([[cluster, m]], index = [cluster], columns = ['sub_cluster', 'subTME'])])

for pathway in gene_df.columns:
    for idx in res.index:
        cur = meta.loc[meta.sub_cluster == idx]
        res.loc[idx, pathway] = gene_df.loc[cur.index, pathway].mean()
res = res.iloc[:, 1:]

res = res.replace('module1', 'subTME-IS')
res = res.replace('module2', 'subTME-PSE')
res = res.replace('module3', 'subTME-ICI')

res.to_csv("/public/home/yuwenqi/sc-data/selected/review/crc/subTME_diff_plot/plot_data.csv")

