setwd("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_all/")
data_all = read.csv("./gsva_all.csv")
data = data_all
bioCol=c("#1c79c0","#ff3399", "#0389ff", "#5529ff", "#ff890f", "#ff000f")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="Proliferation", color = "cluster",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=cluster),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=cluster),position = position_jitterdodge(),size=0.5)
p1=p1+theme(legend.position = 'right', legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave("./gsva_all_major.pdf", p1,device = 'pdf', width = 4.5, height = 3, units = 'in')

data = subset(data_all, cluster == 'Epi')
p1 = ggplot(data,aes(x=ident,y=Proliferation),fill=ident)+
  scale_fill_manual(values=c('#ff890f','#7CFC00','#20B2AA',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#FF1493', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=ident),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=60,hjust = 1.1,vjust=1.05),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "Proliferation")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_all_epi.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')


data = subset(data_all, cluster == 'Mye')
p1 = ggplot(data,aes(x=ident,y=Proliferation),fill=ident)+
  scale_fill_manual(values=c('#ff890f','#7CFC00','#20B2AA',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#FF1493', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=ident),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=60,hjust = 1.1,vjust=1.05),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "Proliferation")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_all_mye.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')

data = subset(data_all, cluster == 'TNK')
p1 = ggplot(data,aes(x=ident,y=Proliferation),fill=ident)+
  scale_fill_manual(values=c('#ff890f','#7CFC00','#20B2AA',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#FF1493', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=ident),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=60,hjust = 1.1,vjust=1.05),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "Proliferation")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_all_tnk.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')

data = subset(data_all, cluster == 'Fib')
p1 = ggplot(data,aes(x=ident,y=Proliferation),fill=ident)+
  scale_fill_manual(values=c('#ff890f','#7CFC00','#20B2AA',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#FF1493', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=ident),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=60,hjust = 1.1,vjust=1.05),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "Proliferation")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_all_fib.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')

data = subset(data_all, cluster == 'End')
p1 = ggplot(data,aes(x=ident,y=Proliferation),fill=ident)+
  scale_fill_manual(values=c('#ff890f','#7CFC00','#20B2AA',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#FF1493', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=ident),show.legend = F,trim = F)+ #outlier.shape = NA
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=60,hjust = 1.1,vjust=1.05),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "Proliferation")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_all_end.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')