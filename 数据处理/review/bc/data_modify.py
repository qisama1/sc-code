import pandas as pd
from sklearn import preprocessing

data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/bc/gene_exp_cluster.csv", index_col = 0)
gene_selected = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/cluster_genes_all_used.csv")
data = data.loc[:, gene_selected.gene]
n = pd.DataFrame(preprocessing.scale(data, axis=0), index = data.index, columns=data.columns)

n.to_csv("/public/home/yuwenqi/sc-data/selected/review/bc/heatmap_plot_data.csv")