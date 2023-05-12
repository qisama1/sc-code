import pandas as pd
corr_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/PDAC/corr_mye.csv", index_col = 0)
corr_caf = pd.read_csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/module/cci_an_sheet/pathway_corr/caf_corr.csv", index_col = 0)
corr_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/PDAC/corr_epi.csv", index_col = 0)
corr_mye_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-score/CRC/corr_gene_gene.csv", index_col = 0)


tam_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/tam_gene.csv")
epi_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/epi_gene.csv")
caf_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/caf_gene.csv")

res_mye = pd.concat([corr_mye.loc[tam_gene.gene, ['M2 Macrophage Polarization']], corr_mye_gene.loc[tam_gene.gene, ['CD163']]], axis=1)
res_epi = corr_epi.loc[epi_gene.gene, ['pEMT']]
res_caf = corr_caf.loc[caf_gene.gene, ['CAF.1']]
res_caf.columns = ['CAF']

res_mye.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/BC/res_mye.csv")
res_epi.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/BC/res_epi.csv")
res_caf.to_csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr_heatmap/BC/res_caf.csv")