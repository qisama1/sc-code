# 按照差异基因选
```python
all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all.csv", index_col = 0)
meta = all[all.cluster == 'Epi']
epi_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Epi/epi_res.csv", index_col=0, encoding = 'gbk')
marker = pd.read_csv("/public/home/yuwenqi/data/Data26/Epi/markers.csv")

sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.sub_cluster))
sg = sg.fillna(0)

sg = sg.replace(' ', '_')

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, epi_res.loc[epi_res.cluster == t.cluster].index] = score
sg = sg.loc[(sg.index.str[0:3] != 'RPL') & (sg.index.str[0:3] != 'RPS') & (sg.index.str[0:3] != 'MT-') & (sg.index.str[0:3] != 'RLS')]
for i in sg.columns:
    t = sg.loc[:, i]
    res = pd.DataFrame(t.loc[t > 0.3].sort_values(ascending = False)[0:30])
    res.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Epi/crc_epi_score/"+ i + "_pct.csv")
```