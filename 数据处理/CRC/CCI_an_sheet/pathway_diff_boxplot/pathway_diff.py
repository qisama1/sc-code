meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/CRC1/all4.csv", index_col=0)
meta_t = meta.loc[~meta['orig.ident'].str.contains('_n')]

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/t.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/epi.csv", index_col = 0).T
gsva_caf = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/caf.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/mye.csv", index_col = 0).T
gsva_end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/gsva/end.csv", index_col = 0).T

gsva_all['cluster'] = meta['cluster']
gsva_all['ident'] = meta['sub_cluster']
gsva_all.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_all/gsva_all.csv")

gsva_t['ident'] = meta['sub_cluster']
gsva_t.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_tnk/gsva_t.csv")

gsva_epi['ident'] = meta['sub_cluster']
gsva_epi.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_epi/gsva_epi.csv")

gsva_caf['ident'] = meta['sub_cluster']
gsva_caf.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_fib/gsva_fib.csv")

gsva_mye['ident'] = meta['sub_cluster']
gsva_mye.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_mye/gsva_mye.csv")


gsva_end['ident'] = meta['sub_cluster']
gsva_end.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_end/gsva_end.csv")