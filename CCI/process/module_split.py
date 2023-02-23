import pandas as pd

v = v.loc[v.sum(axis=1) != 0] 

module1 = pd.DataFrame()
module2 = pd.DataFrame()
module3 = pd.DataFrame()
module4 = pd.DataFrame()
others = pd.DataFrame()

for i in v.columns:
    cluster1 = i.split('|')[0]
    cluster2 = i.split('|')[1]
    if ((cluster1 in set(module.module1.dropna())) & (cluster2 in set(module.module1.dropna()))) :
        module1 = pd.concat([module1, v.loc[:, i]], axis=1)
        continue
    if ((cluster1 in set(module.module2.dropna())) & (cluster2 in set(module.module2.dropna()))) :
        module2 = pd.concat([module2, v.loc[:, i]], axis=1)
        continue
    if ((cluster1 in set(module.module3.dropna())) & (cluster2 in set(module.module3.dropna()))) :
        module3 = pd.concat([module3, v.loc[:, i]], axis=1)
        continue
    others = pd.concat([others, v.loc[:, i]], axis=1)
data = pd.concat([module1, module2, module3], axis=1)