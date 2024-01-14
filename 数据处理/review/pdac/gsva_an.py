import pandas as pd

gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/gsva_score/all2.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/all2.csv", index_col=0)
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/35/module/module.csv", encoding = 'gbk')
res = pd.DataFrame()

for m in ['module1', 'module2', 'module3', 'module4']:
    for cluster in module[m].dropna():
        res = pd.concat([res, pd.DataFrame([[cluster, m]], index = [cluster], columns = ['sub_cluster', 'subTME'])])

for pathway in gsva.columns:
    for idx in res.index:
        cur = meta.loc[meta.ident == idx]
        res.loc[idx, pathway] = gsva.loc[cur.index, pathway].mean()

res = res.replace('module1', 'subTME-PSE')
res = res.replace('module2', 'subTME-PSE')
res = res.replace('module3', 'subTME-IS+MRM')
res = res.replace('module4', 'subTME-MRM')
res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gsva_an/pdac/plot_data.csv")

