import pandas as pd
data = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/concated_data.csv", index_col = 0)
sample = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/concated_sample.csv", index_col = 0)
genes = pd.read_csv("/public/home/yuwenqi/targers_work/data_process/genes.csv")

idx = sample.loc[sample.cancer_type == 'Normal', ].index
data_used = data.loc[idx[idx.isin(data.index)], genes.gene[genes.gene.isin(data.columns)]]
data_used['Tissue'] = sample['Tissue']


cancers = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/check/cancers.csv")


tissue_gene_mean = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/exp_select/tissue_gene_mean.csv", index_col = 0)
tissue_gene_mean = tissue_gene_mean.loc[cancers, ]

for num in range(13,27):
  normal_low_genes = tissue_gene_mean.columns[((tissue_gene_mean <= 5).sum() >= num) & (tissue_gene_mean.loc['LIHC',] <= 5)]
  pd.DataFrame(normal_low_genes, columns=['gene']).to_csv("/public/home/yuwenqi/works/zhangzimei/filter/cancer_" + str(num) +".csv")

data_n_mean = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/val_filter/data_n_mean.csv", index_col = 0)
data_n_mean = data_n_mean.loc[cancers, ]

for cancer in cancers.Tissue:
    t = data_n_mean.loc[cancer, ]
    cur_data = data_n_mean.loc[:,t[t <= 5].index]
    #cur_data.loc[:, cur_data.columns[cur_data.columns.isin(genes.gene)]]
    path = "/public/home/yuwenqi/targers_work/data_concat/val_filter/normal_filter/" + cancer + "/"
    for num in range(13, 26):
        normal_low_genes = cur_data.columns[((cur_data <= 5).sum() >= num)]
        pd.DataFrame(normal_low_genes, columns=['gene']).to_csv(path + str(num) +".csv")

data_t_mean = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/val_filter/data_t_mean.csv", index_col = 0)

data_t_mean = data_t_mean.loc[cancers, ]
for cancer in cancers:
    t = data_t_mean.loc[cancer, ]
    cur_data = data_t_mean.loc[:,t[t >= 5].index]
    path = "/public/home/yuwenqi/targers_work/data_concat/val_filter/tumor_filter/" + cancer + "/"
    for num in range(13, 26):
        normal_low_genes = cur_data.columns[((cur_data >= 5).sum() >= num)]
        pd.DataFrame(normal_low_genes, columns=['gene']).to_csv(path + str(num) +".csv")

for cancer in cancers:
    t = data_t_mean.loc[cancer, ]
    n = data_n_mean.loc[cancer, ]
    cur_data_t = data_t_mean.loc[:,t[t >= 5].index]
    cur_data_n = data_n_mean.loc[:, n[n <= 5].index]
    path = "/public/home/yuwenqi/targers_work/data_concat/val_filter/concat_filter/" + cancer + "/"
    os.mkdir(path)
    for num in range(13, 26):
        tumor_high_genes = set(cur_data_t.columns[(cur_data_t >= 5).sum() >= num])
        normal_low_genes = set(cur_data_n.columns[(cur_data_n <= 5).sum() >= num])
        cur_genes = tumor_high_genes & normal_low_genes
        pd.DataFrame(cur_genes, columns=['gene']).to_csv(path + str(num) +".csv")
