d = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cellchat/module1_df_net.csv", index_col=0)
res = pd.DataFrame(index = set(d.ligand + '|' + d.receptor),columns = set(d.source + "|" + d.target)).fillna(0)


for i in range(len(d)):
    cur = d.iloc[i, ]
    idx = cur.ligand + '|' + cur.receptor
    col = cur.source + '|' + cur.target
    if (cur.pval <= 0.05): res.loc[idx, col] = cur.prob
