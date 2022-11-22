## lc1
mye = pd.read_csv("/public/home/yuwenqi/data/LungData/tLung/cMye4.csv", index_col=0)
num2 = mye.groupby('Sample').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(mye['x']):
  x = mye[mye['x'] == i]
  num1 = x.groupby('Sample').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)
lc1 = t.loc[:, ['Macro_APOE/CTSZ']]
lc1['Type'] = 'LC1'
data = lc1

## lc2
mye = pd.read_csv("/public/home/yuwenqi/data/LungData2/tumor/Mye/submye2.csv", index_col=0)
p = pd.read_csv("/public/home/yuwenqi/data/LungData2/patient.csv", index_col=0)
mye = pd.concat([mye,p.loc[mye.index,]],axis=1)

num2 = mye.groupby('PatientNumber').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(mye['x']):
  x = mye[mye['x'] == i]
  num1 = x.groupby('PatientNumber').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)
lc2 = t.loc[:, ['Macro_APOE/CTSZ']]
lc2['Type'] = 'LC2'
data = pd.concat([lc1, lc2])

## mmrp
d = pd.read_csv("/public/home/yuwenqi/data/Data26/mmrp_tumor_stage_v2.csv", index_col=0)
mye = d.loc[d.ident == 'Mye']
num2 = mye.groupby('pid').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(mye['cell_type']):
  x = mye[mye['cell_type'] == i]
  num1 = x.groupby('pid').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)
mmrp = t.loc[:, ['Macro_APOE/CTSZ']]
mmrp['Type'] = 'MMRp'
data = pd.concat([data, mmrp])

## mmrd
d = pd.read_csv("/public/home/yuwenqi/data/Data26/mmrd_tumor_stage_v2.csv", index_col=0)
mye = d.loc[d.ident == 'Mye']
num2 = mye.groupby('pid').count().iloc[:, 0]
t = pd.DataFrame()
for i in set(mye['cell_type']):
  x = mye[mye['cell_type'] == i]
  num1 = x.groupby('pid').count().iloc[:, 0]
  cur = pd.DataFrame((num1 / num2)).fillna(0)
  cur.columns = [i]
  t = pd.concat([t, cur], axis=1)
mmrd = t.loc[:, ['Macro_APOE/CTSZ']]
mmrd['Type'] = 'MMRd'
data = pd.concat([data, mmrd])