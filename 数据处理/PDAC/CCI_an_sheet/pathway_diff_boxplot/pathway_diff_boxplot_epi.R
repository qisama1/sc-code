setwd("/public/home/yuwenqi/sc-data/selected/35/module/cci_an_sheet/pathway_corr/gsva_diff/gsva_epi/")
data = read.csv("./gsva_epi.csv")
p1 = ggplot(data,aes(x=ident,y=pEMT),fill=ident)+
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
       title = "pEMT")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_epi_pEMT.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')

p1 = ggplot(data,aes(x=ident,y=cEMT),fill=ident)+
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
       title = "cEMT")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_epi_cEMT.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')

p1 = ggplot(data,aes(x=ident,y=Mesenchymal),fill=ident)+
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
       title = "Mesenchymal")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gsva_epi_Mesenchymal.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')



