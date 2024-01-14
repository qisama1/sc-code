```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_used.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all3.csv", row.names=1)

gene = read.csv("/public/home/yuwenqi/sc-data/selected/review/genes_used.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/review/crc/gene_df.csv")
```

```python
import pandas as pd

gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/review/crc/gene_df.csv", index_col = 0)
gene_df = gene_df.T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all3.csv", index_col = 0)

gene_df = gene_df.loc[meta.index, ]
## groupby 大类 mean
gene_df['cluster'] = meta['cluster']

gene_df.groupby('cluster').mean().to_csv("/public/home/yuwenqi/sc-data/selected/review/crc/gene_exp_cluster.csv")
```