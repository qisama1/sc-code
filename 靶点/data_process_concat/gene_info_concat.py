import os
import pandas as pd
dir = "/public/home/yuwenqi/targers_work/data_concat/val_filter/tumor_normal/max_cox/"
file_list = os.listdir(dir)

all_genes = pd.DataFrame()

for file in file_list :
    file_path = dir + "" + file
    genes = pd.read_csv(file_path)
    if len(genes) == 0 :
        continue
    t = pd.DataFrame(genes.gene)
    t.columns = [file.split('.')[0]]
    all_genes = pd.concat([all_genes, t], axis = 1)
