end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/end_meta.csv", index_col=0)
end['cluster'] = 'End'

mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Mye/mye_meta.csv", index_col=0)
mye['cluster'] = 'Mye'

fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Fib/fib_meta.csv", index_col=0)
fib['cluster'] = 'Fib'

tnk = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/TNK/tnk_meta.csv", index_col=0)
tnk['cluster'] = 'TNK'
tnk.columns = tnk.columns.str.replace('sub_clusters', 'sub_cluster')

epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/Epi/epi_meta.csv", index_col=0)
epi['cluster'] = 'Epi'
meta = pd.concat([end, mye, fib, tnk, epi])

## 

# 算所有样本的所有亚类属于其大类的细胞比例
num2 = meta.groupby(['cluster', 'orig.ident']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.sub_cluster):
    x = meta[meta['sub_cluster'] == i]
    num1 = x.groupby('orig.ident').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)
t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/percent.csv")

# 算大类细胞比例
num2 = meta.groupby('orig.ident').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta['cluster']):
  x = meta[meta['cluster'] == i]
  num1 = x.groupby('orig.ident').count().iloc[:, 0]
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
  mm = percent.loc[:, set(meta.loc[meta.cluster == i, ].sub_cluster)]
  mm['type'] = t['type']
  mm.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/percent/" + i + ".csv")
  mm['type1'] = t['type1']
  mm['type2'] = t['type2']
  mm.loc[mm.type == 'tumor'].to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/percent/" + i + "_tumor.csv")

# 算大类数量
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0) 
data=pd.DataFrame(index=list(set(meta['cluster'])), columns=list(set(meta['orig.ident'])))
for i in data.columns:
	for j in Counter(meta[meta['orig.ident']==i]['cluster']):
		data.loc[j, i] = Counter(meta[meta['orig.ident']==i]['cluster'])[j]
type = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/percent/major_type.csv", index_col = 0)
data_t = data.loc[:, type.type == 'tumor']
data_n = data.loc[:, type.type == 'normal']
data_t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/num/major_t.csv")
data_n.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/num/major_n.csv")
for cluster in data.index:
  meta_cur = meta[meta.cluster == cluster]
  data=pd.DataFrame(index=list(set(meta_cur['sub_cluster'])), columns=list(set(meta['orig.ident'])))
  for i in data.columns:
	  for j in Counter(meta_cur[meta_cur['orig.ident']==i]['sub_cluster']):
		  data.loc[j, i] = Counter(meta_cur[meta_cur['orig.ident']==i]['sub_cluster'])[j]
  data_t = data.loc[:, type.type == 'tumor']
  data_n = data.loc[:, type.type == 'normal']    
  data_t.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/num/" + cluster + "_t.csv")
  data_n.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/percent/num/" + cluster + "_n.csv")

