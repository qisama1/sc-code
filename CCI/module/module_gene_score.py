# bc
write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module_genes/genes/module2&3.csv")

write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/genes/module1_module2_genes.csv")
write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module_genes/genes/module2.csv")
write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/genes/module1_genes.csv")

# crc
write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module_genes/genes/module1.csv")

data.loc[(data.cluster1.isin(module.module1)) | (data.cluster2.isin(module.module1))]

genes = pd.DataFrame(list(set(data.gene_a) | set(data.gene_b)))

## 从sce中获取
write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module_genes/genes/module2&3.csv")

write.csv(as.matrix(scRNA[genes$X0, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module_genes/genes/module1.csv")

## res
for i in range(len(data)):
    cur = data.iloc[i, ]
    cluster1 = cur.cluster1
    cluster2 = cur.cluster2
    gene1 = cur.gene_a
    gene2 = cur.gene_b
    genes.loc[tumor.loc[tumor.sub_cluster == cluster1, ].index, [gene1, 'pid']].groupby('pid').sum()
    genes.loc[tumor.loc[tumor.sub_cluster == cluster2, ].index, [gene2, 'pid']].groupby('pid').sum()

sample = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/precent_diff/dataset3/sampleinfo.csv", index_col=0)

tumor_genes = genes.loc[tumor.index,]
tumor_genes['pid'] = tumor.pid
t = pd.DataFrame(index = set(tumor['pid']), columns = set(data.genes))
t = t.fillna(0)
for i in range(len(data)):
    cluster1 = data.iloc[i, ].cluster1
    cluster2 = data.iloc[i, ].cluster2
    geneA = data.iloc[i, ].gene_a
    geneB = data.iloc[i, ].gene_b
    score1 = tumor_genes.loc[tumor.loc[tumor.sub_cluster == cluster1].index, [geneA, 'pid']].groupby(tumor_genes['pid']).sum()
    score2 = tumor_genes.loc[tumor.loc[tumor.sub_cluster == cluster1].index, [geneB, 'pid']].groupby('pid').sum()
    t.loc[score1.index, data.iloc[i, ].genes] += score1.iloc[:, 0].values
    t.loc[score2.index, data.iloc[i, ].genes] += score2.iloc[:, 0].values
t_tumor = t
t_tumor['mean_score'] = t_tumor.mean(axis=1)
t_tumor['type'] = 'Tumor'

normal_genes = genes.loc[normal.index,]
normal_genes['pid'] = normal.pid
t = pd.DataFrame(index = set(normal['pid']), columns = set(data.genes))
t = t.fillna(0)
for i in range(len(data)):
    cluster1 = data.iloc[i, ].cluster1
    cluster2 = data.iloc[i, ].cluster2
    geneA = data.iloc[i, ].gene_a
    geneB = data.iloc[i, ].gene_b
    score1 = normal_genes.loc[normal.loc[normal.sub_cluster == cluster1].index, [geneA, 'pid']].groupby('pid').sum()
    score2 = normal_genes.loc[normal.loc[normal.sub_cluster == cluster1].index, [geneB, 'pid']].groupby('pid').sum()
    t.loc[score1.index, data.iloc[i, ].genes] += score1.iloc[:, 0].values
    t.loc[score2.index, data.iloc[i, ].genes] += score2.iloc[:, 0].values
t_normal = t
t_normal['mean_score'] = t_normal.mean(axis=1)
t_normal['type'] = 'Normal'

pd.concat([t_tumor, t_normal])
t_tumor['grade'] = sample.loc[t_tumor.index, 'grade']