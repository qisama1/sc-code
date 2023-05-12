```py
sub1 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/sub1.csv")
sub2 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/sub2.csv")
sub3 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/sub3.csv")

genes = set()
# 得到需要的基因
for sub in [sub1, sub2, sub3]:
    for i in sub.gene:
        idx = i
        if (len(idx.split('_')) > 1):
            idx = idx.split('_')[0]
            gene_c = i.split('_')[1]
            genes.add(gene_c)
        genes.add(idx.split('|')[0])
        genes.add(idx.split('|')[1])
pd.DataFrame(genes, columns = ['gene'])
# 构建结果列表
cols = []
for sub in [sub1, sub2, sub3]:
    cols.append(sub.type.dropna().values)
idxs = genes

# 构建结果集，应该构建三个
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
res1 = pd.DataFrame(index = genes, columns = sub1.type.dropna())
idxs = set(sub1.gene.dropna().values) | set(sub2.gene.dropna().values) | set(sub3.gene.dropna().values)

res1 = pd.DataFrame(index = idxs, columns = sub1.type.dropna()).fillna(0)
for idx in res1.index:
    for col in res1.columns:
        if (len(idx.split('_')) > 1):
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            gene3 = idx.split("|")[1].split('_')[1]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res1.loc[idx, col] = (data.loc[meta.loc[meta.ident == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene2].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene3].mean()) / 3
        else:
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res1.loc[idx, col] = (data.loc[meta.loc[meta.ident == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene2].mean()) / 2
cellchat = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellchat/module1_df_net_v.csv", index_col = 0)
for i in res1.index:
    for j in res1.columns:
        if (i not in cellchat.index):
            res1.loc[i, j] = 0
            continue
        if (cellchat.loc[i, j] == 0):
            res1.loc[i, j] = 0
res2 = pd.DataFrame(index = idxs, columns = sub2.type.dropna()).fillna(0)
for idx in res2.index:
    for col in res2.columns:
        if (len(idx.split('_')) > 1):
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            gene3 = idx.split("|")[1].split('_')[1]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res2.loc[idx, col] = (data.loc[meta.loc[meta.ident == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene2].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene3].mean()) / 3
        else:
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res2.loc[idx, col] = (data.loc[meta.loc[meta.ident == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene2].mean()) / 2
cellchat = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellchat/module2_df_net_v.csv", index_col = 0)
for i in res2.index:
    for j in res2.columns:
        if (i not in cellchat.index):
            res2.loc[i, j] = 0
            continue
        if (cellchat.loc[i, j] == 0):
            res2.loc[i, j] = 0
res3 = pd.DataFrame(index = idxs, columns = sub3.type.dropna()).fillna(0)
for idx in res3.index:
    for col in res3.columns:
        if (len(idx.split('_')) > 1):
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            gene3 = idx.split("|")[1].split('_')[1]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res3.loc[idx, col] = (data.loc[meta.loc[meta.ident == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene2].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene3].mean()) / 3
        else:
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res3.loc[idx, col] = (data.loc[meta.loc[meta.ident == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.ident == cluster2].index, gene2].mean()) / 2
cellchat = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellchat/module3_df_net_v.csv", index_col = 0)
for i in res3.index:
    for j in res3.columns:
        if (i not in cellchat.index):
            res3.loc[i, j] = 0
            continue
        if (cellchat.loc[i, j] == 0):
            res3.loc[i, j] = 0
sorted_idx = list(res1.loc[(res1.sum(axis=1) > res2.sum(axis=1)) & (res1.sum(axis=1) > res3.sum(axis=1))].index)+ list(res2.loc[(res2.sum(axis=1) > res1.sum(axis=1)) & (res2.sum(axis=1) > res3.sum(axis=1))].index)+ list(res3.loc[(res3.sum(axis=1) > res1.sum(axis=1)) & (res3.sum(axis=1) > res2.sum(axis=1))].index)

sorted_idx.reverse()

pd.concat([res1.loc[sorted_idx, res1.sum().sort_values(ascending=False).index], res2.loc[sorted_idx, res2.sum().sort_values(ascending=False).index], res3.loc[sorted_idx, res3.sum().sort_values(ascending=False).index]], axis=1)
pd.concat([res1.loc[sorted_idx, res1.sum().sort_values(ascending=False).index], res2.loc[sorted_idx, res2.sum().sort_values(ascending=False).index], res3.loc[sorted_idx, res3.sum().sort_values(ascending=False).index]], axis=1).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/bubble_base.csv")

```

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/ESCA/cellchat/gene_df.csv")
```
```