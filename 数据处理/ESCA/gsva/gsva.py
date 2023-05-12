from scipy.stats import spearmanr
gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ssgsea/esca_ssgsea.csv", index_col=0).T

module_score = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/module/tpm_compare/tpm_module_type.csv", index_col = 0)
module_score = module_score.drop('type', axis=1)

corr_res = pd.DataFrame(index = module_score.columns, columns = gsva.columns).fillna(0)

for idx in corr_res.index:
    for col in corr_res.columns:
        s_r = spearmanr(module_score[idx], gsva[col])
        if (s_r[1] <= 0.05):
            corr_res.loc[idx, col] = s_r[0]
    