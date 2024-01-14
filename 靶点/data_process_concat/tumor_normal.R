library(tidyverse)
# 区分TCGA的癌症正常
sample = read.csv("/public/home/yuwenqi/targers_work/data_concat/concated_sample.csv", row.names = 1)

genes = read.csv("/public/home/yuwenqi/targers_work/data_process/genes.csv")

cancers = c('ESCA', 'UCEC', 'HNSC', 'KIRP', 'OV',
'LUSC', 'STAD', 'UCS', 'THCA',  'CESC', 'PRAD',
'COAD', 'SARC', 'LIHC', 'BRCA', 'READ',
'KICH', 'LUAD', 'BLCA', 'PAAD', 'CHOL','THYM',
'ACC', 'GBM', 'SKCM', 'KIRC')

data = qread("/public/home/yuwenqi/targers_work/data_concat/data_used.qs")
data = data[rownames(sample),]

data['cancer_type'] = sample['cancer_type']
data['cancer'] = sample['Tissue']

filtered_data = data[data$cancer %in% cancers, ]
setwd("/public/home/yuwenqi/targers_work/data_concat/pval_select/")
fs = list.files('./', 'csv')


for (x in fs) {
    setwd("/public/home/yuwenqi/targers_work/data_concat/pval_select/")
    genes = read.csv(x)
    c = str_split(x,'[.;]',simplify = T)[,1]
    setwd("/public/home/yuwenqi/targers_work/data_concat/boxplot2/")
    dir.create(c,recursive = T)
    setwd(paste0("./", c))
    for (gene in genes$genes) {
        if (gene %in% colnames(filtered_data)) {
            plot_data = filtered_data[, c(gene, 'cancer_type', 'cancer')]
            colnames(plot_data) = c('gene', 'cancer_type', 'cancer')
            my_boxplot(plot_data, gene)
        }
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
        geom_boxplot(aes(fill=cancer_type),show.legend = T, outlier.colour = NA)+ #outlier.shape = NA
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
        coord_cartesian(ylim = c(0, 6))+
        stat_compare_means(aes(group=cancer_type), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                            label = "p.signif", method='wilcox.test', label.x = 0.55, label.y = 5)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                            
    p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
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


p1 = ggplot(data,aes(x=variable,y=value),fill=cancer_type)+
        scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                    '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                    '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                    '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                    '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                    '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                    '#D9D9D9', '#EEAEEE'))+
        geom_boxplot(aes(fill=cancer_type),show.legend = T, outlier.colour = NA)+ #outlier.shape = NA
        theme(panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(),
                axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
                plot.title = element_text(hjust=0.5),)+
        labs(x=NULL,
            y=NULL,
            title = "")+
        #theme_minimal() 
        #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
        coord_cartesian(ylim = c(0, 50))+
        stat_compare_means(aes(group=cancer_type), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                            label = "p.signif", method='wilcox.test', label.x = 0.55, label.y = 45)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                            
    p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
    ggsave("/public/home/yuwenqi/targers_work/data_concat/ASPM-C/t_n.pdf", width = 6, height = 2.6, units = 'in')