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

data = read.csv("/public/home/yuwenqi/sc-data/selected/review/bc/subTME_diff_plot/plot_data.csv", row.names=1)
data = melt(data, id.vars = c('subTME'))

p1 = ggplot(data,aes(x=variable,y=value),fill=subTME)+
    scale_fill_manual(values=c('#aed4e9','#db6968','#20B2AA',  '#F0F8FF',
                                '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                                '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                                '#AB82FF', '#FF1493', '#B03060', '#708090', 
                                '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                                '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                                '#D9D9D9', '#EEAEEE'))+
    geom_boxplot(aes(fill=subTME),show.legend = T)+ #outlier.shape = NA
    theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(),
            axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
            plot.title = element_text(hjust=0.5),)+
    labs(x=NULL,
        y=NULL,
        title = NULL)+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='wilcox.test', label.x = 0.55, label.y = max(data['value']) - 0.7,size=8)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(1, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 20) )

ggsave("/public/home/yuwenqi/sc-data/selected/review/bc/subTME_diff_plot/subTME_diff_plot.pdf", width = 30, height = 2.6, units = 'in')

