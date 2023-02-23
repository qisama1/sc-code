res = pd.DataFrame(index = set(d.ligand + '|' + d.receptor),columns = set(d.source + "|" + d.target))
res = res.fillna(0)

for i in range(len(d)):
    cur = d.iloc[i, ]
    idx = cur.ligand + '|' + cur.receptor
    col = cur.source + '|' + cur.target
    res.loc[idx, col] = cur.prob

pd.DataFrame(list(set(res.index.str.split('|').str[0]) | set(res.index.str.split('|').str[1])))