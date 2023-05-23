import pandas as pd

gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/all2.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all3.csv", index_col = 0)
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/modules/module.csv", index_col=0)

res = pd.DataFrame()

for m in ['module1', 'module2', 'module3']:
    for cluster in module[m].dropna():
        res = pd.concat([res, pd.DataFrame([[cluster, m]], index = [cluster], columns = ['sub_cluster', 'subTME'])])

for pathway in gsva.columns:
    for idx in res.index:
        cur = meta.loc[meta.sub_cluster == idx]
        res.loc[idx, pathway] = gsva.loc[cur.index, pathway].mean()

res = res.replace('module1', 'subTME-IS')
res = res.replace('module2', 'subTME-PSE')
res = res.replace('module3', 'subTME-ICI')
res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gsva_an/crc/plot_data.csv")

