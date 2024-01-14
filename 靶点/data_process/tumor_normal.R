library(tidyverse)
# 区分TCGA的癌症正常
sample = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/sample.csv", row.names = 1)

['ESCA', 'UCEC', 'HNSC', 'KIRP', 'OV',
'LUSC', 'STAD', 'UCS', 'THCA', 'MESO', 'CESC', 'PRAD',
'COAD', 'SARC', 'LIHC', 'BRCA', 'READ',
'KICH', 'LUAD', 'BLCA', 'PAAD', 'CHOL', 'UVM', 'THYM',
'TGCT', 'ACC', 'GBM', 'SKCM', 'KIRC']

genes = read.csv("/public/home/yuwenqi/targers_work/data_process/genes.csv")

cancers = c('ESCA', 'UCEC', 'HNSC', 'KIRP', 'OV',
'LUSC', 'STAD', 'UCS', 'THCA', 'MESO', 'CESC', 'PRAD',
'COAD', 'SARC', 'LIHC', 'BRCA', 'READ',
'KICH', 'LUAD', 'BLCA', 'PAAD', 'CHOL', 'UVM', 'THYM',
'TGCT', 'ACC', 'GBM', 'SKCM', 'KIRC')


data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/tpm3.csv", row.names=1)
data = data[rownames(sample),]
all_cancers = read.csv("/public/home/yuwenqi/targers_work/data_process/all_cancers.csv")
data['cancer_type'] = factor(ifelse(as.numeric(substr(rownames(data),14,15)) < 10,'Tumor','Normal'))
data['cancer'] = sample['disease_code']
filtered_data = data[data$cancer %in% all_cancers$cancers, ]


for (gene in genes$gene) {
    if (gene %in% colnames(filtered_data)) {
        plot_data = filtered_data[, c(gene, 'cancer_type', 'cancer')]
        colnames(plot_data) = c('gene', 'cancer_type', 'cancer')
        my_boxplot(plot_data, gene)
    }
}

setwd("/public/home/yuwenqi/targers_work/data_process/tumor_normal_plot/")
my_boxplot = function(plot_data, gene_name) {
    p1 = ggplot(plot_data,aes(x=cancer,y=gene),fill=cancer_type)+
        scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                    '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                    '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                    '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                    '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                    '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                    '#D9D9D9', '#EEAEEE'))+
        geom_boxplot(aes(fill=cancer_type),show.legend = T)+ #outlier.shape = NA
        theme(panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(),
                axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
                plot.title = element_text(hjust=0.5),)+
        labs(x=NULL,
            y=NULL,
            title = gene_name)+
        #theme_minimal() 
        #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
        #coord_cartesian(ylim = c(0, 2))+
        stat_compare_means(aes(group=cancer_type), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                            label = "p.signif", method='wilcox.test', label.x = 0.55, label.y = max(plot_data['gene']) - 0.15)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                            
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
    ggsave(paste0("./", gene_name, ".pdf"), width = 14, height = 2.6, units = 'in')
}

#wilcox
setwd("/public/home/yuwenqi/targers_work/data_process/pval_select/")
for (c in all_cancers$cancers) {
    cancer_data = subset(data, cancer == c)
    tumor_data = subset(cancer_data, cancer_type == 'Tumor')
    normal_data = subset(cancer_data, cancer_type == 'Normal')
    gene_selected = data.frame()
    for (gene in genes$gene) {
        if (gene %in% colnames(tumor_data)) {
            tumor_used = tumor_data[, c(gene, 'cancer_type')]
            normal_used = normal_data[, c(gene, 'cancer_type')]
            colnames(tumor_used) = c('gene', 'cancer_type')
            colnames(normal_used) = c('gene', 'cancer_type')
            if (wilcox.test(tumor_used$gene, normal_used$gene)$p.value < 0.05) {
                t = data.frame(gene)
                rownames(t) = c(gene)
                colnames(t) = c('genes')
                gene_selected = rbind(gene_selected, t)
            }
        }
    }
    write.csv(gene_selected, paste0("./", c, ".csv"))
}