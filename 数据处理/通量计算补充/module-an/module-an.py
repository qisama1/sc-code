c1 = pd.read_csv("/public/home/yuwenqi/scFEA/moduleinfo_IN.csv", index_col = 0)
c2 = pd.read_csv("/public/home/yuwenqi/scFEA/moduleinfo_OUT.csv", index_col = 0)

set1 = set(c1['IN'])
set2 = set(c2['OUT'])

set3 = set()
set4 = set()

for i in set1:
  if ('in' not in i) & ('IN' not in i):
	  set3.add(i)
	  
for i in set2:
  if ('out' not in i) & ('OUT' not in i):
	  set4.add(i)

sct = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/module-an/pdac/sct.csv", index_col = 0)

data=[]
for i in list(set3 | set4):
		data.append(- sct.loc[:, c1[c1['IN'] == i].index].sum(axis=1) + sct.loc[:, c2[c2['OUT'] == i].index].sum(axis=1))
		
meta=pd.concat(data, axis=1)
meta.columns=list(set3|set4)