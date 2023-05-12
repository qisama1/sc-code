mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/mye_info.csv", index_col=0)
fib = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/fib_info.csv", index_col=0)
end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/end_info.csv", index_col=0)
b = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/b_info.csv", index_col=0)
tnk = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/tnk_info.csv", index_col=0)
epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/37/epi_info.csv", index_col=0)

mye['cluster'] = 'Mye'

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

# 算大类细胞比例
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all2.csv", index_col=0)
# 每个大类内部的细胞比例
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/percent_new/percent.csv", index_col=0)
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

p_info = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/p_info.csv", index_col=0)
p = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/percent_diff/tumor_normal.csv", index_col=0)
# 每个大类内部的细胞比例

## 切割每个大类

percent.loc[percent.index.str[0] == 'T', 'type'] = 'Tumor'
percent.loc[percent.index.str[0] == 'N', 'type'] = 'Normal'
percent.index = percent.index.str[:-1] 

for i in set(meta.cluster):
  mm = percent.loc[:, set(meta.loc[meta.cluster == i, ].ident)]
  mm['type'] = percent['type']
  mm['stage'] = percent['stage']
  mm.to_csv("/public/home/yuwenqi/sc-data/selected/35/module/percent_stage_compare/" + i + ".csv")
  mm.loc[mm.type == 'Tumor'].to_csv("/public/home/yuwenqi/sc-data/selected/35/module/percent_stage_compare/" + i + "_tumor.csv")


## 切割每个大类
percent = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/percent_new/percent.csv", index_col=0)
percent.loc[percent.index.str[0] == 'T', 'type'] = 'Tumor'
percent.loc[percent.index.str[0] == 'N', 'type'] = 'Normal'
percent.index = percent.index.str[:-1] 

for i in set(meta.cluster):
  mm = percent.loc[:, set(meta.loc[meta.cluster == i, ].sub_cluster)]
  mm['type'] = percent['type']
  mm['grade'] = percent['grade']
  mm['node_type'] = percent['node_type']
  mm.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/percent_new/percent_diff/" + i + ".csv")
  mm.loc[mm.type == 'tumor'].to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/percent_new/percent_diff/" + i + "_tumor.csv")

