data=read.csv("./ppt31/ppt31_mmrp_cxcr6.csv")
p1 = ggplot(data,aes(x=cluster,y=CXCR6),fill=cluster)+
  scale_fill_manual(values=c('#FFDEAD','#7CFC00','#20B2AA',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#FF1493', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=cluster),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "CXCR6")  
#annotate("text", label='***', x = 'Macro_APOE/CTSZ', y = 4)
#theme_minimal() 
#geom_sina(aes(color=cluster), size=0.1,alpha=0.2) 
annotate("text", x = 'CD4Teff', y = 5, label='***', alpha=1)

ggsave(filename = "./ppt31/ppt31_mmrp_CXCR6.pdf", plot = p1, device = 'pdf', width = 2.5, height = 2, units = 'in')

data=read.csv("./ppt31/ppt31_mmrp_cxcl16.csv")
p2 = ggplot(data,aes(x=cluster,y=CXCL16),fill=cluster)+
  scale_fill_manual(values=c('#FFDEAD','#FF1493','#7CFC00',  '#F0F8FF',
                             '#778899', '#EEA2AD', '#E0FFFF', '#CDCD00', 
                             '#C2C2C2', '#C1CDC1', '#8E8E38', '#CD1076',
                             '#AB82FF', '#20B2AA', '#B03060', '#708090', 
                             '#4F94CD', '#473C8B', '#212121', '#00E5EE',
                             '#006400', '#8B0000', '#8B5F65', '#CD8C95',
                             '#D9D9D9', '#EEAEEE'))+
  geom_boxplot(aes(fill=cluster),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "CXCL16")  +
  scale_x_discrete(limits=c('Macro_APOE/CTSZ','Macro_Pro','Macro_CCL19','Macro_CXCL8','Macro_FTL','Mono', 'DC'))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =4.5, label='***', alpha=1)


ggsave(filename = "./ppt31/ppt31_mmrp_CXCL16.pdf", plot = p2, device = 'pdf', width = 2, height = 3, units = 'in')
