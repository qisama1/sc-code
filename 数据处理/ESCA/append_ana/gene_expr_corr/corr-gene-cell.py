data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gene_exp/ESCA/gene_df.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)
data = data.loc[meta.index, ]
data['sample'] = meta['sample']

part1_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/part1_gene.csv")
part1_percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/part1_percent.csv")
part2_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/part2_gene.csv")
part2_percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/part2_percent.csv")
part3_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/part3_gene.csv")
part3_percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/part3_percent.csv")
## percent
num2 = meta.groupby(['cluster', 'sample']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.ident):
    x = meta[meta['ident'] == i]
    num1 = x.groupby('sample').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)

## gene
res = pd.DataFrame(index = t.index, columns = part1_gene.type + "_" + part1_gene.gene).fillna(0)
for i in part1_gene.index:
    v = part1_gene.iloc[i, ]
    g = v.gene
    if (g not in data.columns) :
        continue
    if (v.type == 'TAM'):
        major = meta.loc[meta.loc[meta.ident.str.contains('Macro')].index,]
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'CD8+T'):
        major = meta.loc[meta.loc[meta.ident.str.contains('CD8')].index,]
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d
    if (v.type == 'Epi'):
        major = meta.loc[meta.cluster == 'Epi']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'CAF'):
        major = meta.loc[meta.cluster == 'Fib']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'End'):
        major = meta.loc[meta.cluster == 'End']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]

corr_res = pd.DataFrame(index = res.columns, columns = part1_percent.loc[part1_percent.type == 'ESCA'].cluster)
for i in corr_res.index:
    for j in corr_res.columns:
        if (spearmanr(t[j], res[i])[1] < 0.05):
            corr_res.loc[i, j] = spearmanr(t[j], res[i])[0]
corr_res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/ESCA/part1_corr.csv")

res = pd.DataFrame(index = t.index, columns = part2_gene.type + "_" + part2_gene.gene).fillna(0)
for i in part2_gene.index:
    v = part2_gene.iloc[i, ]
    g = v.gene
    if (g not in data.columns) :
        continue
    if (v.type == 'TAM'):
        major = meta.loc[meta.loc[meta.ident.str.contains('Macro')].index,]
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'CD8+T'):
        major = meta.loc[meta.loc[meta.ident.str.contains('CD8')].index,]
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d
    if (v.type == 'Epi'):
        major = meta.loc[meta.cluster == 'Epi']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'CAF'):
        major = meta.loc[meta.cluster == 'Fib']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'End'):
        major = meta.loc[meta.cluster == 'End']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]

corr_res = pd.DataFrame(index = res.columns, columns = part2_percent.loc[part2_percent.type == 'ESCA'].cluster)
for i in corr_res.index:
    for j in corr_res.columns:
        if (spearmanr(t[j], res[i])[1] < 0.05):
            corr_res.loc[i, j] = spearmanr(t[j], res[i])[0]
corr_res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/ESCA/part2_corr.csv")

res = pd.DataFrame(index = t.index, columns = part3_gene.type + "_" + part3_gene.gene).fillna(0)
for i in part3_gene.index:
    v = part3_gene.iloc[i, ]
    g = v.gene
    if (g not in data.columns) :
        continue
    if (v.type == 'TAM'):
        major = meta.loc[meta.loc[meta.ident.str.contains('Macro')].index,]
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'CD8+T'):
        major = meta.loc[meta.loc[meta.ident.str.contains('CD8')].index,]
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d
    if (v.type == 'Epi'):
        major = meta.loc[meta.cluster == 'Epi']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'CAF'):
        major = meta.loc[meta.cluster == 'Fib']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]
    if (v.type == 'End'):
        major = meta.loc[meta.cluster == 'End']
        d = data.loc[major.index, ].groupby('sample').mean().loc[:, [g]]
        res.loc[d.index, v.type + "_" + v.gene] = d[g]

corr_res = pd.DataFrame(index = res.columns, columns = part3_percent.loc[part3_percent.type == 'ESCA'].cluster)
for i in corr_res.index:
    for j in corr_res.columns:
        if (spearmanr(t[j], res[i])[1] < 0.05):
            corr_res.loc[i, j] = spearmanr(t[j], res[i])[0]
corr_res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/ESCA/part3_corr.csv")