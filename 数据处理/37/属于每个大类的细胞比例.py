mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/mye_info.csv", index_col=0)
fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/fib_info.csv", index_col=0)
tnk = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/tnk_info.csv", index_col=0)
epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/epi_info.csv", index_col=0)

mye['cluster'] = 'Mye'
fib['cluster'] = 'Fib'
tnk['cluster'] = 'TNK'
epi['cluster'] = 'Epi'

# 算所有样本的所有亚类属于其大类的细胞比例
num2 = meta.groupby(['cluster', 'sample']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.ident):
    x = meta[meta['ident'] == i]
    num1 = x.groupby('sample').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)

# 算大类细胞比例
num2 = meta.groupby('sample').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta['cluster']):
  x = meta[meta['cluster'] == i]
  num1 = x.groupby('sample').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)
t = t.fillna(0)
# 每个大类内部的细胞比例

## 切割每个大类

percent.loc[percent.index.str[0] == 'T', 'type'] = 'Tumor'
percent.loc[percent.index.str[0] == 'N', 'type'] = 'Normal'
percent.index = percent.index.str[:-1] 

for i in set(meta.cluster):
  mm = percent.loc[:, set(meta.loc[meta.cluster == i, ].ident)]
  mm['Grade'] = percent['Grade']
  mm['Stage'] = percent['Stage']
  mm['stage_node'] = percent['stage_node']
  mm.to_csv("/public/home/yuwenqi/sc-data/selected/37/module/percent_stage_compare/" + i + ".csv")

