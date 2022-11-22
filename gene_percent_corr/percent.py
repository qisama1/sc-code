import pandas as pd

expr = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/gene_percent_corr/genes.csv", index_col = 0)
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/bc3_meta_t_major2.csv", index_col = 0)

# gene_expr_mean
expr['cluster'] = meta['orig.ident']
expr_mean_patient = expr.groupby('cluster').mean()

# percent
num2 = meta.groupby('orig.ident').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta['sub_clusters']):
  x = meta[meta['sub_clusters'] == i]
  num1 = x.groupby('orig.ident').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)

percent_patient = t
