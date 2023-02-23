percent_patient
expr_mean_patient

from scipy.stats import spearmanr

## ecm_myCAF:Cycling
gene_name = ['COL1A1', 'CXCL12', 'FN1', 'THBS1', 'THBS2', 'GRN', 'GRN', 'CXCL12', 'CXCL14']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Cycling"]))

gene_name = ['SDC1', 'CXCR4', 'SDC1', 'SDC1', 'SDC1', 'SORT1', 'TNFRSF1A', 'CXCR4', 'CXCR4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["ecm_myCAF"]))

## ecm_myCAF:ML_1
gene_name = ['COL1A1', 'CXCL12', 'FN1', 'THBS1', 'THBS2', 'GRN', 'GRN', 'CXCL12', 'CXCL14']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["ML_1"]))

gene_name = ['SDC1', 'CXCR4', 'SDC1', 'SDC1', 'SDC1', 'SORT1', 'TNFRSF1A', 'CXCR4', 'CXCR4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["ecm_myCAF"]))


## Angiogenic_EC:Cycling
gene_name = ['JAM2', 'HSPG2', 'CXCL12', 'FN1', 'COL4A1', 'GRN', 'JAG1', 'GRN', 'GRN', 'CXCL12']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Cycling"]))

gene_name = ['F11R', 'DAG1', 'CXCR4', 'SDC1', 'SDC1', 'SORT1', 'CD46', 'TNFRSF1A', 'TNFRSF1B', 'CXCR4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Angiogenic_EC"]))

## DC_1|ProT
spearmanr(expr_mean_patient["TIGIT"], percent_patient["DC_1"])

## Angiogenic_EC:ProT
gene_name = ['ICAM2', 'ICAM1', 'COL4A1', 'NECTIN2', 'NECTIN2', 'PECAM1']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["ProT"]))

gene_name = ['ITGAL', 'ITGAL', 'CD44', 'TIGIT', 'TIGIT', 'CD38']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Angiogenic_EC"]))

## DC_1:Treg
gene_name = ['CD86', 'TNF', 'CXCL9', 'TNFSF13B', 'TNFSF13B', 'TNF', 'TNF', 'P2RY6', 'CD86', 'CD86']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Treg"]))

gene_name = ['CTLA4', 'VSIR', 'CXCR3', 'HLA-DPB1', 'TFRC', 'ICOS', 'TNFRSF1B', 'NAMPT', 'CD28', 'CTLA4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["DC_1"]))

## Cycling:ecm_myCAF
spearmanr(expr_mean_patient["MDK"], percent_patient["ecm_myCAF"])

gene_name = ['SDC1', 'SDC2', 'LRP1']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Cycling"]))

gene_name = ['SDC1', 'SDC2', 'LRP1']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["ML_1"]))

## ecm_myCAF:Macro_SPP1
gene_name = ['LAMB2', 'COL1A1', 'FN1', 'FN1', 'CXCL12', 'CXCL12', 'CXCL14']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Macro_SPP1"]))

gene_name = ['ITGA2', 'ITGA2', 'ITGAV', 'ITGA5', 'DPP4', 'CXCR4', 'CXCR4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["ecm_myCAF"]))

## Mono_1:Angiogenic_EC
gene_name = ['ITGB2', 'TNFSF13B', 'TNFSF13B']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Angiogenic_EC"]))

gene_name = ['ICAM1', 'HLA-DPB1', 'CD40']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Mono_1"]))

## Cycling|Angiogenic_EC
gene_name = ['EFNA1', 'EFNA3', 'EFNA4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Angiogenic_EC"]))

gene_name = ['EPHA4', 'EPHA2']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Cycling"]))    

## Angiogenic_EC:CD8Teff
gene_name = ['NECTIN2', 'SELP', 'SELE']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["CD8Teff"]))

gene_name = ['TIGIT', 'SELPLG', 'SELPLG']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Angiogenic_EC"]))

## Angiogenic_EC:Mono_1
gene_name = ['JAG1', 'PDGFB', 'CCL14', 'CXCL12']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Mono_1"]))

gene_name = ['CD46', 'LRP1', 'CCR1', 'CXCR4']
for i in gene_name:
    print(spearmanr(expr_mean_patient[i], percent_patient["Angiogenic_EC"]))