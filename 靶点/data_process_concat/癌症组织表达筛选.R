library(qs)
library(ggsci)
library(ggplot2)
library(cowplot)
library(data.table)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(tidyverse)

data = qread("/public/home/yuwenqi/targers_work/tpm_data.qs")

all_cancers = read.csv("/public/home/yuwenqi/targers_work/data_process/all_cancers.csv")
sample = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/sample.csv", row.names = 1)
data['cancer_type'] = factor(ifelse(as.numeric(substr(rownames(data),14,15)) < 10,'Tumor','Normal'))
data['cancer'] = sample['disease_code']


setwd("/public/home/yuwenqi/targers_work/data_process/pval_select/")
fs = list.files('./', 'csv')
selected_normal = read.csv("/public/home/yuwenqi/targers_work/data_process/selected_data/cancers_gene_selected.csv")
all_cancers_plot_data = data.frame()

for (x in fs) {
    c = str_split(x,'[.;]',simplify = T)[,1]
    cancer_data = subset(data, cancer == c)
    setwd("/public/home/yuwenqi/targers_work/data_process/pval_select/")
    genes = read.csv(x)
    
    tumor_data = subset(cancer_data, cancer_type == 'Tumor')
    normal_data = subset(cancer_data, cancer_type == 'Normal')
    tumor_data = tumor_data[, genes$genes]
    tumor_data = tumor_data[, selected_normal[selected_normal[, c] != "", c]]
    
    temp_data = tumor_data[,]
    tumor_data['cancer_type'] = cancer_data[rownames(tumor_data), 'cancer_type']
    plot_data = melt(tumor_data, id.vars = c('cancer_type'))
    plot_data['cancer'] = c
    setwd("/public/home/yuwenqi/targers_work/data_process/exp_plot_tumor/")
    all_cancers_plot_data = rbind(all_cancers_plot_data, plot_data)
    exp_boxplot(plot_data, c)
    selected_data = temp_data[,colMeans(temp_data) > quantile(plot_data$value, 0.6)]
    setwd("/public/home/yuwenqi/targers_work/data_process/selected_data_tumor/")
    write.csv(selected_data, paste0("./", c, ".csv"))
    print(dim(selected_data))
}


exp_boxplot = function(plot_data, cancer) {
    p1 = ggplot(plot_data,aes(x=cancer_type,y=value),fill=cancer_type)+
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
            title = cancer)                           
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
    ggsave(paste0("./", cancer, ".pdf"), width = 1.2, height = 2.6, units = 'in')
}


all_exp_boxplot = function(plot_data) {
    p1 = ggplot(plot_data,aes(x=cancer,y=value),fill=cancer)+
        scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                    '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                    '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                    '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                    '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                    '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                    '#D9D9D9', '#EEAEEE'))+
        geom_boxplot(aes(fill=cancer),show.legend = T)+ #outlier.shape = NA
        theme(panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(),
                axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
                plot.title = element_text(hjust=0.5),)+
        labs(x=NULL,
            y=NULL,
            title = "all_cancers")                           
    p1=p1+theme(legend.position = "none", legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
    ggsave(paste0("./", "all", ".pdf"), width = 6, height = 2.6, units = 'in')
}
setwd("/public/home/yuwenqi/targers_work/data_process/exp_plot_tumor/")
all_exp_boxplot(all_cancers_plot_data)