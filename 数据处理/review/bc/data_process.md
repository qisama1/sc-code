```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", row.names=1)
gene = read.csv("/public/home/yuwenqi/sc-data/selected/review/genes_used.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/review/bc/gene_df.csv")
```

```python
import pandas as pd

gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/bc/gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col = 0)

gene_df = gene_df.loc[meta.index, ]
## groupby 大类 mean
gene_df['cluster'] = meta['cluster']

gene_df.groupby('cluster').mean().to_csv("/public/home/yuwenqi/sc-data/selected/review/bc/gene_exp_cluster.csv")
```