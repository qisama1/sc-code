## bc
num2 = meta.groupby(['major', 'orig.ident']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(d2['cluster']):
  x = d2[d2['cluster'] == i]
  num1 = x.groupby('pid').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)

num2 = meta.groupby(['major', 'orig.ident']).count().iloc[:, 0]
module_all = pd.DataFrame()
for m in module.columns:
    t = pd.DataFrame()
    for i in set(module[m].dropna()):
        x = meta[meta['sub_clusters'] == i]
        if len(x) == 0 :
            continue
        num1 = x.groupby('orig.ident').count().iloc[:, 0]
        cur = pd.DataFrame((num1 / (num2.loc[x.major[0]])))
        cur.columns = [i]
        t = pd.concat([t, cur], axis=1)
    t = t.fillna(0)
    module_v = pd.DataFrame(t.sum(axis = 1) / len(t.columns))
    module_v.columns = [m]
    module_all = pd.concat([module_all, module_v], axis=1)

sample = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/precent_diff/dataset3/sampleinfo.csv", index_col=0)

num2 = normal.groupby(['major', 'orig.ident']).count().iloc[:, 0]
module_all = pd.DataFrame()
for m in module.columns:
    t = pd.DataFrame()
    for i in set(module[m].dropna()):
        x = normal[normal['sub_clusters'] == i]
        if len(x) == 0 :
            continue
        num1 = x.groupby('orig.ident').count().iloc[:, 0]
        cur = pd.DataFrame((num1 / (num2.loc[x.major[0]])))
        cur.columns = [i]
        t = pd.concat([t, cur], axis=1)
    t = t.fillna(0)
    module_v = pd.DataFrame(t.sum(axis = 1) / len(t.columns))
    module_v.columns = [m]
    module_all = pd.concat([module_all, module_v], axis=1)

## crc
tumor = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_tumor.csv", index_col=0)

num2 = tumor.groupby(['cluster', 'pid']).count().iloc[:, 0]
module_tumor = pd.DataFrame()
for m in module.columns:
    t = pd.DataFrame()
    for i in set(module[m].dropna()):
        x = tumor[tumor['sub_cluster'] == i]
        if len(x) == 0 :
            continue
        num1 = x.groupby('pid').count().iloc[:, 0]
        cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
        cur.columns = [i]
        t = pd.concat([t, cur], axis=1)
    t = t.fillna(0)
    module_v = pd.DataFrame(t.sum(axis = 1) / len(t.columns))
    module_v.columns = [m]
    module_tumor = pd.concat([module_tumor, module_v], axis=1)

normal = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/crc1_normal.csv", index_col=0)
num2 = normal.groupby(['cluster', 'pid']).count().iloc[:, 0]
module_normal = pd.DataFrame()
for m in module.columns:
    t = pd.DataFrame()
    for i in set(module[m].dropna()):
        x = normal[normal['sub_cluster'] == i]
        if len(x) == 0 :
            continue
        num1 = x.groupby('pid').count().iloc[:, 0]
        cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
        cur.columns = [i]
        t = pd.concat([t, cur], axis=1)
    t = t.fillna(0)
    module_v = pd.DataFrame(t.sum(axis = 1) / len(t.columns))
    module_v.columns = [m]
    module_normal = pd.concat([module_normal, module_v], axis=1)