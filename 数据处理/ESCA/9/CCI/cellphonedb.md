# 处理p-v
```python
p = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/out/pvalues.txt", sep='\t', index_col=0)
v = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/out/means.txt", sep = '\t', index_col=0)

v.index = v.gene_a + '|' + v.gene_b
v = v.loc[~v.index.isna()].iloc[:, 10:]
p.index = p.gene_a + '|' + p.gene_b
p = p.loc[~p.index.isna()].iloc[:, 10:]
p = p[~p.index.duplicated()]
v = v[~v.index.duplicated()]

set(v.columns.str.split('|').str[0])

for idx in v.index:
    for col in v.columns:
        if p.loc[idx, col] > 0.05:
            v.loc[idx, col] = 0

## cnt
cnt = pd.DataFrame(index = set(v.columns.str.split('|').str[0]), columns = set(v.columns.str.split('|').str[0])).fillna(0)

for col in v.columns:
    t = v.loc[:, col]
    cluster1 = col.split('|')[0]
    cluster2 = col.split('|')[1]
    cnt.loc[cluster1, cluster2] = t.loc[t > 0].count()
p.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/p.csv")
v.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/v.csv")
cnt.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/cnt.csv")
```

```R
library(pheatmap)
setwd("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/")
count_matrix = read.csv("/public/home/yuwenqi/sc-data/selected/9/module/cellphonedb/module3/cnt.csv", row.names=1)

show_rownames = T
show_colnames = T
scale="none"
cluster_cols = T
border_color='white'
cluster_rows = T
fontsize_row=11
fontsize_col = 11
main = ''
treeheight_row=0
family='Arial'
treeheight_col = 0
col1 = "dodgerblue4"
col2 = 'peachpuff'
col3 = 'deeppink4'
meta_sep='\t'
pvalues_sep='\t'
pvalue=0.05
col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )

h = pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col)
ggsave(filename = "./cnt_module3_plot.pdf", plot = h, device = 'pdf', width = 6, height = 6, units = 'in')


```