```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
genes = c('GLUL')
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT11-dotplot/bc/gene_df.csv")
```

```python
import pandas as pd
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
meta_t = meta.loc[meta.Tissue == 'Tumor']
major = meta_t.loc[meta_t.sub_cluster.str.contains('Macro')].index

df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT10-dotplot/bc/gene_df.csv", index_col = 0).T
df = df.loc[major, ['GLUL', 'GSS', 'GGT6', 'GGT1', 'GGT7', 'GGT5']]

sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/bc/sct2.csv", index_col = 0)
sct = sct.loc[major, ['M_48', 'M_25', 'M_26']]

metabolism = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/bc/metabolism_t2.csv", index_col = 0)
metabolism = metabolism.loc[major, ['Glutamine']]

gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/all2.csv", index_col = 0).T
gsva = gsva.loc[major, ['M2 Macrophage Polarization', 
'Anti-inflammatory in myeloid cells']]
gsva.columns = ['M2', 'Anti.inflammatory']
pd.concat([df, sct, metabolism, gsva], axis = 1).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT10-dotplot/bc/plot_data.csv")
```
