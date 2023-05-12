```py
sub1 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/sub1.csv")
sub2 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/sub2.csv")
sub4 = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/sub4.csv")
genes = set()
# 得到需要的基因
for sub in [sub1, sub2, sub4]:
    for i in sub.gene.dropna():
        idx = i
        if (len(idx.split('_')) > 1):
            idx = idx.split('_')[0]
            gene_c = i.split('_')[1]
            genes.add(gene_c)
        genes.add(idx.split('|')[0])
        genes.add(idx.split('|')[1])
pd.DataFrame(genes, columns = ['gene'])
pd.DataFrame(genes, columns = ['gene']).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/gene_needed.csv")
# 构建结果列表
cols = []
for sub in [sub1, sub2, sub3]:
    cols.append(sub.type.dropna().values)
idxs = genes

# 构建结果集，应该构建三个
data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
idxs = set(sub1.gene.dropna().values) | set(sub2.gene.dropna().values) | set(sub4.gene.dropna().values)

res1 = pd.DataFrame(index = idxs, columns = sub1.type.dropna()).fillna(0)
for idx in res1.index:
    for col in res1.columns:
        if (len(idx.split('_')) > 1):
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            gene3 = idx.split("|")[1].split('_')[1]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res1.loc[idx, col] = (data.loc[meta.loc[meta.sub_cluster == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene2].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene3].mean()) / 3
        else:
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res1.loc[idx, col] = (data.loc[meta.loc[meta.sub_cluster == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene2].mean()) / 2
cellchat = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cellchat/cellchat_module1_v.csv", index_col = 0)

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
            res2.loc[idx, col] = (data.loc[meta.loc[meta.sub_cluster == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene2].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene3].mean()) / 3
        else:
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res2.loc[idx, col] = (data.loc[meta.loc[meta.sub_cluster == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene2].mean()) / 2
cellchat = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cellchat/cellchat_module2_v.csv", index_col = 0)

for i in res2.index:
    for j in res2.columns:
        if (i not in cellchat.index):
            res2.loc[i, j] = 0
            continue
        if (cellchat.loc[i, j] == 0):
            res2.loc[i, j] = 0

res4 = pd.DataFrame(index = idxs, columns = sub4.type.dropna()).fillna(0)
for idx in res4.index:
    for col in res4.columns:
        if (len(idx.split('_')) > 1):
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            gene3 = idx.split("|")[1].split('_')[1]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res4.loc[idx, col] = (data.loc[meta.loc[meta.sub_cluster == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene2].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene3].mean()) / 3
        else:
            gene1 = idx.split('|')[0]
            gene2 = idx.split("|")[1].split('_')[0]
            cluster1 = col.split('|')[0]
            cluster2 = col.split("|")[1]
            res4.loc[idx, col] = (data.loc[meta.loc[meta.sub_cluster == cluster1].index, gene1].mean() + data.loc[meta.loc[meta.sub_cluster == cluster2].index, gene2].mean()) / 2
cellchat = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cellchat/cellchat_module4_v.csv", index_col = 0)

for i in res4.index:
    for j in res4.columns:
        if (i not in cellchat.index):
            res4.loc[i, j] = 0
            continue
        if (cellchat.loc[i, j] == 0):
            res4.loc[i, j] = 0

sorted_idx = list(res1.loc[(res1.sum(axis=1) > res2.sum(axis=1)) & (res1.sum(axis=1) > res4.sum(axis=1))].index)+ list(res2.loc[(res2.sum(axis=1) > res1.sum(axis=1)) & (res2.sum(axis=1) > res4.sum(axis=1))].index)+ list(res4.loc[(res4.sum(axis=1) > res1.sum(axis=1)) & (res4.sum(axis=1) > res2.sum(axis=1))].index)

sorted_idx.reverse()

pd.concat([res1.loc[sorted_idx, res1.sum().sort_values(ascending=False).index], res2.loc[sorted_idx, res2.sum().sort_values(ascending=False).index], res4.loc[sorted_idx, res4.sum().sort_values(ascending=False).index]], axis=1).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/bubble_base.csv")

```

```R
scRNA = qread("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/bc3.qs")
gene = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/gene_needed.csv")
write.csv(as.matrix(scRNA[gene$gene, ]@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellchat/gene_df.csv")
```
```