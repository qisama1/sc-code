```py
map_info = pd.read_csv("/public/home/yuwenqi/targers_work/getx/map_info.csv", index_col = 0)
map_info = map_info.fillna('nope')

sample = pd.read_csv("/public/home/yuwenqi/targers_work/getx/annotations_v8_GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", index_col = 0, sep='\t')


map_data = pd.DataFrame()
map_sample = pd.DataFrame()
for tcga_name in map_info.index:
    gtex_name = map_info.loc[tcga_name, 'GTEX']
    if (gtex_name == 'nope'):
        continue
    idxs = gtex_sample.loc[gtex_sample.SMTS == gtex_name].index
    cur_matrix = gtex_matrix.loc[idxs, ]
    cur_matrix.index = cur_matrix.index + "_" + tcga_name
    cur_sample = gtex_sample.loc[idxs, ]
    cur_sample.index = cur_sample.index + "_" + tcga_name
    map_data = pd.concat([map_data, cur_matrix])
    map_sample = pd.concat([map_sample, cur_sample])

map_sample['Tissue'] = map_sample.index.str.split('_').str[1]

```

```py
gtex_sample = pd.read_csv("/public/home/yuwenqi/targers_work/getx/annotations_v8_GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", index_col = 0, sep='\t')
tcga_sample = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/sample.csv", index_col = 0)
gtex_matrix = matrix.loc[:, matrix.columns.isin(gtex_sample.index)]
gtex_matrix = gtex_matrix.T

c1 = matrix.columns[matrix.columns.str.contains('TCGA')]
tcga_matrix = matrix.loc[:, c1]
tcga_matrix = tcga_matrix.T
tcga_matrix.index = tcga_matrix.index + 'A'
tcga_matrix = tcga_matrix.loc[tcga_matrix.index.isin(tcga_sample.index), :]

gtex_sample = gtex_sample.loc[gtex_matrix.index, ]
SMTS
```