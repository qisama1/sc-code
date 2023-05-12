meta = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", index_col=0)

gsva_all = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/all.csv", index_col = 0).T
gsva_t = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/t.csv", index_col = 0).T
gsva_epi = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/epi.csv", index_col = 0).T
gsva_caf = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/caf.csv", index_col = 0).T
gsva_mye = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/mye.csv", index_col = 0).T
gsva_end = pd.read_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/gsva_score/end.csv", index_col = 0).T

gsva_all['cluster'] = meta['cluster']
gsva_all['ident'] = meta['ident']
gsva_all.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_all/gsva_all.csv")

gsva_t['ident'] = meta['ident']
gsva_t.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_t/gsva_t.csv")

gsva_epi['ident'] = meta['ident']
gsva_epi.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_epi/gsva_epi.csv")

gsva_caf['ident'] = meta['ident']
gsva_caf.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_fib/gsva_fib.csv")

gsva_mye['ident'] = meta['ident']
gsva_mye.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_mye/gsva_mye.csv")


gsva_end['ident'] = meta['ident']
gsva_end.to_csv("/public/home/yuwenqi/sc-data/selected/9/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_end/gsva_end.csv")