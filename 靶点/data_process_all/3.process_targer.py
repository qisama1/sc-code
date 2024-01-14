import os
import pandas as pd
dir = "/public/home/yuwenqi/targers_work/data_concat/pval_select/"
file_list = os.listdir(dir)

all_genes = pd.DataFrame()

for file in file_list :
    file_path = dir + "" + file
    genes = pd.read_csv(file_path)
    if len(genes) == 0 :
        continue
    t = pd.DataFrame(genes.genes)
    t.columns = [file.split('.')[0]]
    all_genes = pd.concat([all_genes, t], axis = 1)

cancers = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/check/cancers.csv")

data = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/concated_data.csv", index_col = 0)
sample = pd.read_csv("/public/home/yuwenqi/targers_work/data_concat/concated_sample.csv", index_col = 0)
sample = sample.loc[data.index, ]
path = "/public/home/yuwenqi/targers_work/data_concat/val_filter/tumor_normal/normal/"
genes = pd.read_csv("/public/home/yuwenqi/targers_work/data_process/genes.csv")

for cancer in cancers.Tissue:
    cur_path = path + cancer + "/25.csv"
    cur_data = pd.read_csv(cur_path)
    len(cur_data.genes[cur_data.genes.isin(genes.gene)])
    selected_data_n = data.loc[sample[(sample.Tissue == cancer) & (sample.cancer_type =='Normal')].index, ].loc[:, cur_data.genes]
    selected_data = data.loc[sample[(sample.Tissue == cancer) & (sample.cancer_type =='Tumor')].index, ].loc[:, cur_data.genes]
    selected_data.loc[:, 'TIGIT'].mean()
    selected_data_n.loc[:, 'TIGIT'].mean()
    selected_genes = selected_data.loc[:, selected_data.mean() > 5].columns
    t = pd.DataFrame(selected_genes[selected_genes.isin(genes.gene)], columns = ['gene'])
    save_dir = "/public/home/yuwenqi/targers_work/data_concat/val_filter/tumor_normal/max_cox/"
    t.to_csv(save_dir + cancer + ".csv")

    
 
   