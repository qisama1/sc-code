# 表达量筛选
```py
import pandas as pd
data = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/concated_data.csv", index_col = 0)
sample = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/concated_sample.csv", index_col = 0)
genes = pd.read_csv("/public/home/yuwenqi/targers_work/data_process/genes.csv")

idx = sample.loc[sample.cancer_type == 'Normal', ].index
data_used = data.loc[idx[idx.isin(data.index)], genes.gene[genes.gene.isin(data.columns)]]
data_used['Tissue'] = sample['Tissue']

cancers = ['ACC', 'BLCA', 'BRCA', 'COAD', 'ESCA', 'GBM', 'HNSC', 'KICH',
       'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD',
       'PRAD', 'READ', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC',
       'UCS']

tissue_gene_mean = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/exp_select/tissue_gene_mean.csv", index_col = 0)

tissue_gene_mean = tissue_gene_mean.loc[cancers,:]

normal_low_genes = tissue_gene_mean.columns[(tissue_gene_mean <= 5).sum() == len(tissue_gene_mean.index)]

pd.DataFrame(normal_low_genes, columns=['gene']).to_csv("/public/home/yuwenqi/targers_work/data_concat/exp_select/result_genes.csv")
```

# p值筛选

```R
data = qread("/public/home/yuwenqi/targers_work/data_concat/data_used.qs")
sample = read.csv("/public/home/yuwenqi/targers_work/data_concat/concated_sample.csv", row.names = 1)
data[is.na(data)] = 0   
data['cancer_type'] = sample[rownames(data), 'cancer_type']
data['cancer'] = sample[rownames(data), 'Tissue']
cancers = c('ACC', 'BLCA', 'BRCA', 'COAD', 'ESCA', 'GBM', 'HNSC', 'KICH',
       'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD',
       'PRAD', 'READ', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC',
       'UCS')
genes = read.csv("/public/home/yuwenqi/targers_work/data_concat/exp_select/result_genes.csv")
setwd("/public/home/yuwenqi/targers_work/data_concat/pval_select/")
for (c in cancers) {
    cancer_data = subset(data, cancer == c)
    tumor_data = subset(cancer_data, cancer_type == 'Tumor')
    normal_data = subset(cancer_data, cancer_type == 'Normal')
    gene_selected = data.frame()
    for (gene in genes$gene) {
        if (gene %in% colnames(tumor_data)) {
            tumor_used = tumor_data[, c(gene, 'cancer_type')]
            normal_used = normal_data[, c(gene, 'cancer_type')]
            colnames(tumor_used) = c('gene', 'cancer_type')
            colnames(normal_used) = c('gene', 'cancer_type')
            if (!is.na(wilcox.test(tumor_used$gene, normal_used$gene)$p.value) && wilcox.test(tumor_used$gene, normal_used$gene)$p.value < 0.05 && mean(tumor_used[, 'gene']) > 2 * mean(normal_used[, 'gene'])) {
                t = data.frame(gene)
                rownames(t) = c(gene)
                colnames(t) = c('genes')
                gene_selected = rbind(gene_selected, t)
            }
        }
    }
    print(dim(gene_selected))
    write.csv(gene_selected, paste0("./", c, ".csv"))
}
```