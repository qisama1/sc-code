clin_list = pd.read_csv("/public/home/yuwenqi/targers_work/data_process/clin_list.txt", sep='\t', index_col = 0)

tcga_data = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/split_data/UCS_data.csv", index_col = 0)


gene_list = []
for gene in clin_list.index:
    idx = gene.index('(')
    name = gene[idx + 1:-1]
    if name not in tcga_data.columns:
        gene_list.append(gene)

market_list = pd.read_csv("/public/home/yuwenqi/targers_work/data_process/market_list.txt", sep='\t', index_col = 0)



for gene in market_list.index:
    idx = gene.index('(')
    name = gene[idx + 1:-1]
    if name not in tcga_data.columns:
        gene_list.append(gene)

pd.DataFrame(gene_list, columns = ['gene']).to_csv("/public/home/yuwenqi/targers_work/data_process/protein_s.csv")

# 挑选实体瘤
{'ESCA', 'UCEC', 'HNSC', 'KIRP', 'PCPG', 'OV',
'LUSC', 'STAD', 'UCS', 'THCA', 'MESO', 'CESC', 'PRAD',
'COAD', 'SARC', 'LGG', 'LAML', 'LIHC', 'BRCA', 'READ',
'KICH', 'LUAD', 'BLCA', 'PAAD', 'CHOL', 'UVM', 'THYM',
'TGCT', 'ACC', 'GBM', 'SKCM', 'KIRC'}

['ESCA', 'UCEC', 'HNSC', 'KIRP', 'OV',
'LUSC', 'STAD', 'UCS', 'THCA', 'MESO', 'CESC', 'PRAD',
'COAD', 'SARC', 'LIHC', 'BRCA', 'READ',
'KICH', 'LUAD', 'BLCA', 'PAAD', 'CHOL', 'UVM', 'THYM',
'TGCT', 'ACC', 'GBM', 'SKCM', 'KIRC']