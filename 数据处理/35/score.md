# 取基因表达
```R
end = qread("/public/home/yuwenqi/sc-data/selected/35/end.qs")
tnk = qread("/public/home/yuwenqi/sc-data/selected/35/tnk.qs")
epi = qread("/public/home/yuwenqi/sc-data/selected/35/epi.qs")
b = qread("/public/home/yuwenqi/sc-data/selected/35/bcell.qs")
mye = qread("/public/home/yuwenqi/sc-data/selected/35/mye.qs")
fib = qread("/public/home/yuwenqi/sc-data/selected/35/fib.qs")
scRNA = merge(end, tnk)
scRNA = merge(scRNA, epi)
scRNA = merge(scRNA, b)
scRNA = merge(scRNA, mye)
scRNA = merge(scRNA, fib)
scRNA = qread("...")
write.csv(as.matrix(scRNA@assays$RNA@data), "/public/home/yuwenqi/sc-data/selected/35/module/cluster_sg_scores/gene_df.csv")
```

# 算分数
```python
meta = pd.read_csv("all.csv")

s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.ident)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.ident == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/mye/9_mye_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/tnk/9_tnk_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'B'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/bcell/9_b_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/fib/9_fib_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/epi/9_epi_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/endo/9_end_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cluster_sg_scores/s_celltype.csv")
```