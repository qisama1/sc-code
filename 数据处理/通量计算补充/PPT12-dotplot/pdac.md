```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/35/all3.qs")
genes = c('GLUL')
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT12-dotplot/pdac/gene_df.csv")
```

```python
import pandas as pd
df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT12-dotplot/pdac/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]
num2 = meta_t.groupby(['cluster', 'sample']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta_t.ident):
    x = meta_t[meta_t['ident'] == i]
    num1 = x.groupby('sample').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)

major = meta_t.loc[meta_t.ident.str.contains('Macro'), ].index
df = df.loc[major, ]
df.loc[:, ['sample']] = meta_t.loc[df.index, ['sample']]
pd.concat([t['CD8Tex'], df.groupby('sample').mean().fillna(0)], axis = 1).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT12-dotplot/pdac/plot_data.csv")
```

```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT12-dotplot/pdac/plot_data.csv")
sp1 = ggplot(data, aes(x = CD8Tex, y = GLUL))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = T, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) 
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT12-dotplot/pdac/plot.pdf", plot = sp1, device = 'pdf', width = 2.5, height = 2, units = 'in')
```