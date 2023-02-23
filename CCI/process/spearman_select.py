from scipy.preprocessing import spearmanr

num2 = meta.groupby(['major', 'orig.ident']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.sub_clusters):
    x = meta[meta['sub_clusters'] == i]
    if len(x) == 0 :
        continue
    num1 = x.groupby('orig.ident').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.major[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)

from scipy.stats import spearmanr

for i in v.index:
    for j in v.columns:
        gene = i.split('|')[0]
        cluster = j.split('|')[1]
        if (cluster not in set(meta.sub_clusters)):
            v.loc[i, j] = 0
            continue
        if (gene not in data.columns):
            v.loc[i, j] = 0
            continue
        percent = t[cluster]
        gene_score = data.groupby('pid').mean()[gene]
        spearmanr(percent, gene_score)[1]
        if (spearmanr(percent, gene_score)[1] > 0.05) :
            v.loc[i, j] = 0

v.to_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cellchat/select_v.csv")
## crc
num2 = meta.groupby(['cluster', 'pid']).count().iloc[:, 0]
t = pd.DataFrame()
for i in set(meta.sub_cluster):
    x = meta[meta['sub_cluster'] == i]
    num1 = x.groupby('pid').count().iloc[:, 0]
    cur = pd.DataFrame((num1 / (num2.loc[x.cluster[0]])))
    cur.columns = [i]
    t = pd.concat([t, cur], axis=1)
t = t.fillna(0)

for i in v.index:
    for j in v.columns:
        gene = i.split('|')[0]
        cluster = j.split('|')[1]
        if (cluster not in set(meta.sub_cluster)):
            v.loc[i, j] = 0
            continue
        if (gene not in data.columns):
            v.loc[i, j] = 0
            continue
        percent = t[cluster]
        gene_score = data.groupby('pid').mean()[gene]
        if (spearmanr(percent, gene_score)[1] > 0.05) :
            v.loc[i, j] = 0
