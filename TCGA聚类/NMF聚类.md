## 预处理提取tcga数据
```python
```

```R
tcga = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/coad_tcga.csv", row.names=1)
gene = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/needed_gene.csv")
tcga_use = tcga[, gene$Gene]


normalize <- function(x){return((x-min(x))/(max(x)-min(x)))}

EMT_expr2 <- normalize(tcga_use)

setwd("/public/home/yuwenqi/targers_work/data_concat/val_filter/cluster/")
library(NMF)
estimate <- nmf(EMT_expr2, rank=2:10, method="brunet", nrun=100, seed=123)#运行时间久
pdf('NMF-rank.pdf',width = 8,height = 6)
plot(estimate)
dev.off() 

pdf('consensusmap.pdf',width = 15,height = 15,onefile = F)
consensusmap(estimate,annRow = NA,annCol = NA,main = "Consensus matrix",info =FALSE)
dev.off()

rank <- 2
seed <- 123
estimate2 <- nmf(EMT_expr2, rank=2, method="brunet", nrun=100, seed=seed)

group <- predict(estimate2)
group <- as.data.frame(group)
group$group <- paste0('Cluster',group$group)
group$sample <- rownames(group)
new_group = group[order(group$group),]
table(new_group$group)


pre_heatdata <- EMT_expr2[,new_group$sample]

#annotation_col

group2 <- as.data.frame(new_group[,1])
rownames(group2) <- rownames(new_group)
colnames(group2) <- 'group'

#color
annColors <- list()
annColors[['group']] <- c('Cluster1'='red','Cluster2'='steelblue')
#plot
library(pheatmap)
#加深颜色对比度
pre_heatdata <- as.data.frame(t(scale(t(pre_heatdata))))
pre_heatdata[pre_heatdata > 1] <- 1
pre_heatdata[pre_heatdata < -1] <- -1
pdf('EMT_pheatmap.pdf',width = 10,height = 10)
pheatmap(
    pre_heatdata,
    color = colorRampPalette(c('navy','white','firebrick3'))(1000),
    annotation_col = group2,annotation_colors = annColors,treeheight_row = 30,#gaps_col = 126,
    show_rownames = F, show_colnames = F,cluster_rows = T,cluster_cols = F
)
dev.off()

# 后续再求os
```