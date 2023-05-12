genes = set()
an_pair = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellphonedb/an_pair.csv")
for i in an_pair.index:
    v = an_pair.loc[i, ]
    idx = v.gene
    if (len(idx.split('_')) > 1):
        gene_c = idx.split('_')[1]
        idx = idx.split('_')[0]
        genes.add(gene_c)
    genes.add(idx.split('|')[0])
    genes.add(idx.split('|')[1])
pd.DataFrame(genes, columns = ['gene']).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellphonedb/gene_needed.csv")

genes = set()
an_pair = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellchat/an_pair.csv")
for i in an_pair.index:
    v = an_pair.loc[i, ]
    idx = v.gene
    if (len(idx.split('_')) > 1):
        gene_c = idx.split('_')[1]
        idx = idx.split('_')[0]
        genes.add(gene_c)
    genes.add(idx.split('|')[0])
    genes.add(idx.split('|')[1])
pd.DataFrame(genes, columns = ['gene']).to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/CRC/cellchat/gene_needed.csv")