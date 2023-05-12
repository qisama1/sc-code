# 按照差异基因选
```python
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/end_meta.csv", index_col = 0)
end_res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/end_res.csv", index_col=0, encoding = 'gbk')
marker = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/end_markers_batched3.csv")
end_res.index = end_res.index.str.replace('-', '_').str.replace(' ', '_')
end_res.index = end_res.index.str.replace('γ', 'y')


sg = pd.DataFrame(index = set(marker.gene), columns = set(meta.sub_cluster))
sg = sg.fillna(0)
sg = sg.loc[(sg.index.str[0:3] != 'RPL') & (sg.index.str[0:3] != 'RPS') & (sg.index.str[0:3] != 'MT-') & (sg.index.str[0:3] != 'RLS')]
sg.columns = sg.columns.str.replace(' ', '_').str.replace('-', '_')

for i in marker.index:
    t = marker.loc[i,]
    score = 1 - t['pct.2'] / t['pct.1']
    sg.loc[t.gene, end_res.loc[end_res.cluster == t.cluster].index] = score

for i in sg.columns:
    t = sg.loc[:, i]
    res = pd.DataFrame(t.loc[t > 0.3].sort_values(ascending = False)[0:30])
    res.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/End/bc_end_score/"+ i + "_pct.csv")
```