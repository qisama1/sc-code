import pandas as pd
gene_df = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/gene_df.csv", index_col = 0).T
an_pair = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/an_pair.csv")
meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/bc3_meta_all2.csv", index_col=0)

res = pd.DataFrame(index = an_pair.gene, columns = ['ligand_T', 'ligand_N','receptor_T', 'receptor_N'])

for idx in an_pair.index:
    v = an_pair.loc[idx,]
    gene = v.gene
    type = v.type
    cluster1 = type.split("|")[0]
    cluster2 = type.split("|")[1]
    flag = False
    if (len(gene.split('_')) > 1):
        gene_c = gene.split('_')[1]
        gene = gene.split('_')[0]
        flag = True
    gene_a = gene.split('|')[0]
    gene_b = gene.split('|')[1]
    if (cluster1 == 'CAF'):
        idx1 = meta.loc[meta.cluster == 'Fib'].index
        idx1_T = meta.loc[idx1,].loc[meta.Tissue == 'Tumor'].index
        idx1_N = meta.loc[idx1,].loc[meta.Tissue == 'Normal'].index
    elif (cluster1 == 'Epi'):
        idx1 = meta.loc[meta.cluster == 'Epi'].index
        idx1_T = meta.loc[idx1,].loc[meta.Tissue == 'Tumor'].index
        idx1_N = meta.loc[idx1,].loc[meta.Tissue == 'Normal'].index
    elif (cluster1 == 'Macro'):
        idx1 = meta.loc[meta.sub_cluster.str.contains('Macro')].index
        idx1_T = meta.loc[idx1,].loc[meta.Tissue == 'Tumor'].index
        idx1_N = meta.loc[idx1,].loc[meta.Tissue == 'Normal'].index
    elif (cluster1 == 'CD8Teff'):
        idx1 = meta.loc[meta.sub_cluster.str.contains("CD8")].index
        idx1_T = meta.loc[idx1,].loc[meta.Tissue == 'Tumor'].index
        idx1_N = meta.loc[idx1,].loc[meta.Tissue == 'Normal'].index
    else:
        idx1 = meta.loc[meta.sub_cluster == cluster1].index
        idx2_T = meta.loc[idx2,].loc[meta.Tissue == 'Tumor'].index
        idx2_N = meta.loc[idx2,].loc[meta.Tissue == 'Normal'].index
    if (cluster2 == 'CAF'):
        idx2 = meta.loc[meta.cluster == 'Fib'].index
        idx2_T = meta.loc[idx2,].loc[meta.Tissue == 'Tumor'].index
        idx2_N = meta.loc[idx2,].loc[meta.Tissue == 'Normal'].index
    elif (cluster2 == 'Epi'):
        idx2 = meta.loc[meta.cluster == 'Epi'].index
        idx2_T = meta.loc[idx2,].loc[meta.Tissue == 'Tumor'].index
        idx2_N = meta.loc[idx2,].loc[meta.Tissue == 'Normal'].index
    elif (cluster2 == 'Macro'):
        idx2 = meta.loc[meta.sub_cluster.str.contains('Macro')].index
        idx2_T = meta.loc[idx2,].loc[meta.Tissue == 'Tumor'].index
        idx2_N = meta.loc[idx2,].loc[meta.Tissue == 'Normal'].index
    elif (cluster2 == 'CD8Teff'):
        idx2 = meta.loc[meta.sub_cluster.str.contains("CD8")].index
        idx2_T = meta.loc[idx2,].loc[meta.Tissue == 'Tumor'].index
        idx2_N = meta.loc[idx2,].loc[meta.Tissue == 'Normal'].index
    else:
        idx2 = meta.loc[meta.sub_cluster == cluster2].index
        idx2_T = meta.loc[idx2,].loc[meta.Tissue == 'Tumor'].index
        idx2_N = meta.loc[idx2,].loc[meta.Tissue == 'Normal'].index
    res.loc[v.gene, 'ligand_T'] = gene_df.loc[idx1_T, gene_a].mean()
    res.loc[v.gene, 'ligand_N'] = gene_df.loc[idx1_N, gene_a].mean()
    if (flag):
        res.loc[v.gene, 'receptor_T'] = (gene_df.loc[idx2_T, gene_b].mean() + gene_df.loc[idx2_T, gene_c].mean()) / 2
        res.loc[v.gene, 'receptor_N'] = (gene_df.loc[idx2_N, gene_b].mean() + gene_df.loc[idx2_N, gene_c].mean()) / 2
    else :
        res.loc[v.gene, 'receptor_T'] = (gene_df.loc[idx2_T, gene_b].mean())
        res.loc[v.gene, 'receptor_N'] = (gene_df.loc[idx2_N, gene_b].mean())

res.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/DEG_cellphonedb.csv") 