
```py
tcga = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/tcga_tpm_log2.csv", index_col = 0).T
gsva = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/ssgsea/pdac_ssgsea.csv", index_col = 0).T
gsva.index = gsva.index.str.replace('.', '-')
res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/tpm_module_type.csv", index_col = 0)

part1 = gsva.loc[:, ['Matrix remodeling', 'EMT', 'Immune Suppression by Myeloid Cells']]

an_gene = pd.read_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/tcga_heatmap_all/an_gene.csv")

part2 = pd.DataFrame(index = gsva.index, columns = an_gene.gene)
for i in an_gene.index:
    v = an_gene.loc[i, ]
    idx = v.gene
    if (len(idx.split('_')) > 1):
        gene_c = idx.split('_')[1]
        idx = idx.split('_')[0]
        gene_a = idx.split('|')[0]
        gene_b = idx.split('|')[1]
        part2[v] = pd.DataFrame((tcga.loc[gsva.index, gene_a] + tcga.loc[gsva.index, gene_b] + tcga.loc[gsva.index, gene_c]) / 3, columns = [v])
    else: 
        gene_a = idx.split('|')[0]
        gene_b = idx.split('|')[1]
        part2[v] = pd.DataFrame((tcga.loc[gsva.index, gene_a] + tcga.loc[gsva.index, gene_b]) / 2, columns = [v])

data = pd.concat([part1, part2], axis=1)
plot_data = pd.concat([res['subTME3'], data], axis=1)

plot_data.columns = ['subTME-IS', 'Matrix remodeling', 'EMT',
       'Immune Suppression by Myeloid Cells', 'ANXA1|FPR3', 'C3|C3AR1',
       'C5AR1|RPS19', 'CD47|SIRPG', 'CSF1|SIRPA', 'LGALS9|HAVCR2', 'MDK|LRP1',
       'NR3C1|CXCL8', 'PLAU|PLAUR', 'SIRPA|CD47', 'THBS2|CD47',
       'THY1|ITGAM_ITGB2', 'THY1|ITGAX_ITGB2', 'TYROBP|CD44',
       'COL1A1|ITGA2_ITGB1', 'COL1A1|ITGA3_ITGB1', 'COL1A2|ITGA2_ITGB1',
       'COL1A2|ITGA3_ITGB1', 'COL4A2|ITGA2_ITGB1', 'COL4A2|ITGA2_ITGB1',
       'COL6A2|ITGA2_ITGB1', 'COL6A3|ITGA2_ITGB1', 'LAMA4|ITGA2_ITGB1',
       'LAMA4|ITGA2_ITGB1', 'LAMB1|ITGA2_ITGB1', 'LAMC1|ITGA2_ITGB1',
       'THBS2|CD47', 'CCL4|SLC7A1', 'CCL3|IDE']
plot_data = plot_data.loc[res.type == 'Tumor']
plot_res = plot_data.loc[plot_data.sort_values('subTME-IS',ascending=False).index,].T
import sklearn.preprocessing as preprocessing
plot_res = pd.DataFrame(preprocessing.scale(plot_res, axis=1), index = plot_res.index, columns=plot_res.columns)
plot_res[~plot_res.index.duplicated()].to_csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/tcga_heatmap_all/tcga-plot_data.csv")
```

```R
data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/tcga_heatmap_all/tcga-plot_data.csv", row.names=1)
# an = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/tcga_heatmap_all/an.csv", row.names = 1)
# colnames(data) = rownames(an)

bk <- c(seq(min(data), -0.001, (0 - min(data)) / 100),seq(0, max(data), (max(data) - 0.001) / 100))

heatmap = pheatmap(data, 
                   #annotation_row = annotation, # ?????????????????????????????????
                   #annotation_col=an, # ?????????????????????????????????
                   show_colnames = FALSE, # ??????????????????
                   #show_rownames=TRUE,  # ??????????????????
                   fontsize= 10, # ????????????
                   breaks = bk,
                   color = c(colorRampPalette(c("#0c00e6"))(length(bk) / 4), colorRampPalette(c("#0c00e6","#FFFFFF"))(length(bk) / 4), colorRampPalette(colors=c('#FFFFFF', '#e00000'))(length(bk) / 4), colorRampPalette(colors=c('#e00000'))(length(bk) / 4)), # ?????????????????????
                   #color = colorRampPalette(c("#FFFFFF","#e00000"))(100),
                   annotation_legend=TRUE, # ??????????????????
                   border_color=NA,  # ???????????? NA????????????
                   scale="none",  # ???????????????????????????"row"??????????????????"column"??????????????????"none"?????????
                   cluster_rows = FALSE, # ??????????????????
                   cluster_cols = FALSE, # ??????????????????
                   legend_breaks = c(min(data), 0, max(data) - 0.01),
                   legend_labels = c(round(min(data), 1),0,round(max(data) - 0.01, 1)),
                   cellheight = 15,
                   cellwidth = 1,
                   angle_col = 45,
                   fontsize_col = 10,
                   fontsize_row = 10,
                   #gaps_row = c(5,9, 15, 19, 21, 24, 25)
                   #gaps_row = 1:length(rownames(data)),
                   #gaps_col = c(4, 8, 12, 16)
                   
)
pdf("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/tcga_heatmap_all/tcga-plot_data.pdf", width = 10, height = 10)
heatmap
dev.off()
```
