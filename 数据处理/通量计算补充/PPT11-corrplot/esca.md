```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
genes = c('MAF','PLTP','CD163','MS4A4A','CCL18','MSR1','PSAP','APOE','GPNMB','CTSZ','LGMN','CTSL','CTSB', 'CTSD', 'GLUL', 'NPL')
write.csv(as.matrix(scRNA[genes, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT11-corrplot/esca/gene_df.csv")
```

```python
import pandas as pd
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
meta_t = meta.loc[meta['sample'].str.contains('T')]
major = meta_t.loc[meta_t.ident.str.contains('Macro')].index
df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT11-corrplot/esca/gene_df.csv", index_col = 0).T
def spearmanr_pval(x,y):
    if (spearmanr(x, y)[1] > 0.05):
        return 0
    else :
        return spearmanr(x, y)[0]
df.loc[major, ].corr(method = spearmanr_pval).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT11-corrplot/esca/plot_data.csv")
```

```R
library(corrplot)
library(tidyverse)
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT11-corrplot/esca/plot_data.csv", row.names=1)
data = as.matrix(data)

pdf("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT11-corrplot/esca/plot.pdf", width=5, height=5)
col = c(colorRampPalette(c("#0c00e6","#FFFFFF"))(100), colorRampPalette(colors=c('#EEA9B8', '#ff0c00'))(100))
p1=corrplot(data, type = "lower", order = "hclust", tl.col = "black", tl.srt = 45, col = col)
dev.off()
```