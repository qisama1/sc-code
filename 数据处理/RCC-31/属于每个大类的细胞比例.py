mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/mye_info.csv", index_col=0)
fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/fib_info.csv", index_col=0)
end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/end_info.csv", index_col=0)
b = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/b_info.csv", index_col=0)
tnk = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/tnk_info.csv", index_col=0)
epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/19/epi_info.csv", index_col=0)

mye['cluster'] = 'Mye'
fib['cluster'] = 'Fib'
end['cluster'] = 'End'
b['cluster'] = 'B'
tnk['cluster'] = 'tnk'
epi['cluster'] = 'Epi'

meta[meta.Sample_Origin.isin(['nLung', 'tLung'])]

meta = pd.concat([mye, fib, end, b, tnk, epi])
meta_use = meta[meta.Sample_Origin.isin(['nLung', 'tLung'])]
meta = meta_use
# 算所有样本的所有亚类属于其大类的细胞比例
num2 = meta.groupby([ 'cluster','Sample']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.ident):
    x = meta[meta['ident'] == i]
    num1 = x.groupby('Sample').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)

# 算大类细胞比例
num2 = meta.groupby('Sample').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta['cluster']):
  x = meta[meta['cluster'] == i]
  num1 = x.groupby('Sample').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)
t = t.fillna(0)
# 每个大类内部的细胞比例

## 切割每个大类
percent.loc[percent.index.str.contains('T'), 'type'] = 'Tumor'
percent.loc[~percent.index.str.contains('T'), 'type'] = 'Normal'

for i in set(meta.cluster):
  mm = percent.loc[:, set(meta.loc[meta.cluster == i, ].ident)]
  mm['type'] = percent['type']
  mm.to_csv("/public/home/yuwenqi/sc-data/selected/19/module/percent_stage_compare/" + i + ".csv")

for i in set(meta.cluster):
  mm = percent.loc[:, set(meta.loc[meta.cluster == i, ].ident)]
  mm['type'] = percent['type']
  mm['stage'] = percent['stage']
  mm.to_csv("/public/home/yuwenqi/sc-data/selected/35/module/percent_stage_compare/" + i + ".csv")
  mm.loc[mm.type == 'Tumor'].to_csv("/public/home/yuwenqi/sc-data/selected/35/module/percent_stage_compare/" + i + "_tumor.csv")

