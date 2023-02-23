genes = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/genes/module2_genes.csv", index_col=0)
genes = genes.T
genes['pid'] = meta['orig.ident']

t = pd.DataFrame(index = set(meta['orig.ident']), columns = set(module.genes))
t = t.fillna(0)
for i in range(len(module)):
    cluster1 = module.iloc[i, ].cluster1
    cluster2 = module.iloc[i, ].cluster2
    geneA = module.iloc[i, ].gene_a
    geneB = module.iloc[i, ].gene_b
    score1 = genes.loc[meta.loc[meta.sub_clusters == cluster1].index, [geneA, 'pid']].groupby(genes['pid']).sum()
    score2 = genes.loc[meta.loc[meta.sub_clusters == cluster1].index, [geneB, 'pid']].groupby('pid').sum()
    t.loc[score1.index, module.iloc[i, ].genes] += score1.iloc[:, 0].values
    t.loc[score2.index, module.iloc[i, ].genes] += score2.iloc[:, 0].values

t = pd.DataFrame(index = set(meta['orig.ident']), columns = set(module.genes))
t = t.fillna(0)
for i in range(len(module)):
    cluster1 = module.iloc[i, ].cluster1
    cluster2 = module.iloc[i, ].cluster2
    geneA = module.iloc[i, ].gene_a
    geneB = module.iloc[i, ].gene_b
    score1 = genes.loc[meta.loc[meta.sub_clusters == cluster1].index, [geneA, 'pid']].groupby(genes['pid']).sum()
    score2 = genes.loc[meta.loc[meta.sub_clusters == cluster1].index, [geneB, 'pid']].groupby('pid').sum()
    t.loc[score1.index, module.iloc[i, ].genes] += score1.iloc[:, 0].values
    t.loc[score2.index, module.iloc[i, ].genes] += score2.iloc[:, 0].values



num2 = d2.groupby('pid').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(d2['cluster']):
  x = d2[d2['cluster'] == i]
  num1 = x.groupby('pid').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)