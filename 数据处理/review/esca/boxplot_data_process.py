import pandas as pd
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/esca/gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col = 0)
meta.loc[meta['sample'].str.contains("N"), 'Tissue'] = 'Normal'
meta.loc[~meta['sample'].str.contains("N"), 'Tissue'] = 'Tumor'

genes = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/cluster_genes.csv")
gene_df = gene_df.loc[meta.index,:]
clusters = ['Mye', 'TNK', 'End', 'Fib', 'Epi']

for cluster in clusters:
    t = gene_df.loc[:, genes.loc[:, cluster].dropna()[genes.loc[:, cluster].dropna().isin(gene_df.columns)]]
    t['type'] = meta['Tissue']
    t = t.loc[meta[meta.cluster == cluster].index]
    t.to_csv("/public/home/yuwenqi/sc-data/selected/review/esca/tumor_normal_boxplot/" + cluster + ".csv")   



