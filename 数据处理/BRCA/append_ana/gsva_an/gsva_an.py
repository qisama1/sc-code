import pandas as pd

gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/gsva/all2.csv", index_col = 0).T
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)
module = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc_modules.csv", encoding = 'gbk')

pd.DataFrame([['Macro_APOE', 'module']], columns = ['sub_cluster', 'subTME'])

res = pd.DataFrame()

for m in ['module1', 'module2', 'module3', 'module4']:
    for cluster in module[m].dropna():
        res = pd.concat([res, pd.DataFrame([[cluster, m]], index = [cluster], columns = ['sub_cluster', 'subTME'])])

for pathway in gsva.columns:
    for idx in res.index:
        cur = meta.loc[meta.sub_cluster == idx]
        res.loc[idx, pathway] = gsva.loc[cur.index, pathway].mean()


res = res.replace('module1', 'subTME-MRM')
res = res.replace('module2', 'subTME-IS')
res = res.replace('module3', 'intermediate')
res = res.replace('module4', 'subTME-ICI')

res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/gsva_an/bc/plot_data.csv")

