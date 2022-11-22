import pandas as pd

# cellphonedb
v = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cellphonedb/res/v2.csv", index_col=0)
p = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cellphonedb/res/p2.csv", index_col=0)


ligand = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/nichenet/res/CD8Tem_1/CD8Tem_1_potential.csv",index_col=0)

res = v.loc[:, (v.columns.str.contains('Cycling')) & ((v.columns.str.contains('Angiogenic_EC')))]

res['gene_a'] = res.index.str.split('|').str[0]
res['gene_b'] = res.index.str.split('|').str[1]

res_lr = res.loc[(res.gene_a.isin(ligand.columns)) | (res.gene_b.isin(ligand.columns)),]
res_lr = res_lr.loc[res_lr.sum(axis=1) > 0.2, ]

res_lr.to_csv()

# cellchat

cellchat_v = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cellchat/df_v2.csv", index_col = 0)
cellchat_p = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/cellchat/df_p2.csv", index_col = 0)

cellchat_v = cellchat_v.T
#cellchat_v.loc[:,cellchat_v.columns.str.contains('Mono_1')]

sender_rec = ['Angiogenic_EC:Cycling']

res = cellchat_v.loc[:, sender_rec]

res['gene_a'] = res.index.str.split('_').str[0]
res['gene_b'] = res.index.str.split('_').str[1]

res_lr = res.loc[(res.gene_a.isin(ligand.columns)) | (res.gene_b.isin(ligand.columns)),]
res_lr = res_lr.loc[res_lr.sum(axis=1) > 0.01, ]

