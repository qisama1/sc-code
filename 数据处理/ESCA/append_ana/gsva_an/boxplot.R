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

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/gsva_an/crc/plot_data.csv", row.names=1)
setwd("/public/home/yuwenqi/sc-data/selected/append_ana/gsva_an/crc/")

p1 = ggplot(data,aes(x=subTME,y=Checkpoint.molecules),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Checkpoint.molecules')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Checkpoint.molecules']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Checkpoint.molecules.pdf", width = 1.3, height = 2.3, units = 'in')

p1 = ggplot(data,aes(x=subTME,y=Myeloid.cells.traffic),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Myeloid.cells.traffic')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Myeloid.cells.traffic']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Myeloid.cells.traffic.pdf", width = 1.3, height = 2.3, units = 'in')

p1 = ggplot(data,aes(x=subTME,y=Macrophage.and.DC.traffic),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Macrophage.and.DC.traffic')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Macrophage.and.DC.traffic']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Macrophage.and.DC.traffic.pdf", width = 1.3, height = 2.3, units = 'in')

 
p1 = ggplot(data,aes(x=subTME,y=Matrix.remodeling),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Matrix.remodeling')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Matrix.remodeling']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Matrix.remodeling.pdf", width = 1.3, height = 2.3, units = 'in') 
 
p1 = ggplot(data,aes(x=subTME,y=Angiogenesis),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Angiogenesis')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Angiogenesis']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Angiogenesis.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=Tumor.proliferation),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Tumor.proliferation')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Tumor.proliferation']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Tumor.proliferation.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=Epithelial.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Epithelial.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Epithelial.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Epithelial.cells.pdf", width = 1.3, height = 2.3, units = 'in')    

p1 = ggplot(data,aes(x=subTME,y=T.cell.Exhaustion),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'T.cell.Exhaustion')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['T.cell.Exhaustion']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./T.cell.Exhaustion.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=Interferon.endothelial.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Interferon.endothelial.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Interferon.endothelial.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Interferon.endothelial.cells.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=Angiogenic.endothelial.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Angiogenic.endothelial.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Angiogenic.endothelial.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Angiogenic.endothelial.cells.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=pEMT),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'pEMT')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['pEMT']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./pEMT.pdf", width = 1.3, height = 2.3, units = 'in') 
        
p1 = ggplot(data,aes(x=subTME,y=myCAF),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'pEMT')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['myCAF']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./myCAF.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=Anti.inflammatory.in.T.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'pEMT')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Anti.inflammatory.in.T.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Anti.inflammatory.in.T.cells.pdf", width = 1.3, height = 2.3, units = 'in') 

p1 = ggplot(data,aes(x=subTME,y=Anergy.in.T.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Anergy.in.T.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Anergy.in.T.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Anergy.in.T.cells.pdf", width = 1.3, height = 2.3, units = 'in') 

p1 = ggplot(data,aes(x=subTME,y=Pro.inflammatory.in.T.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Pro.inflammatory.in.T.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Pro.inflammatory.in.T.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Pro.inflammatory.in.T.cells.pdf", width = 1.3, height = 2.3, units = 'in') 

p1 = ggplot(data,aes(x=subTME,y=M1.Macrophage.Polarization),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'M1.Macrophage.Polarization')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['M1.Macrophage.Polarization']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./M1.Macrophage.Polarization.pdf", width = 1.3, height = 2.3, units = 'in') 


p1 = ggplot(data,aes(x=subTME,y=M2.Macrophage.Polarization),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'M2.Macrophage.Polarization')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['M2.Macrophage.Polarization']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./M2.Macrophage.Polarization.pdf", width = 1.3, height = 2.3, units = 'in')

p1 = ggplot(data,aes(x=subTME,y=Pro.inflammatory.in.myeloid.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Pro.inflammatory.in.myeloid.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Pro.inflammatory.in.myeloid.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Pro.inflammatory.in.myeloid.cells.pdf", width = 1.3, height = 2.3, units = 'in')

p1 = ggplot(data,aes(x=subTME,y=Anti.inflammatory.in.myeloid.cells),fill=subTME)+
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
            plot.title = element_text(hjust=0.5),) +
    labs(x=NULL,
        y=NULL,
        title = 'Anti.inflammatory.in.myeloid.cells')+
    #theme_minimal() 
    #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
    #coord_cartesian(ylim = c(0, 2))+
    stat_compare_means(aes(group=subTME), symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", method='kruskal.test', label.x = 0.55, label.y = max(data['Anti.inflammatory.in.myeloid.cells']) - 0.009)  #+ stat_summary(aes(group=type), fun='mean', geom="point")
                        
p1=p1+theme(legend.background = element_blank(), legend.key = element_blank() ,legend.key.size = unit(0.2, "cm"), legend.title = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA), text = element_text(size = 6) )
p1 = p1 + guides(fill=FALSE)
ggsave("./Anti.inflammatory.in.myeloid.cells.pdf", width = 1.3, height = 2.3, units = 'in')