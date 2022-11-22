# 散点图
data = read.csv("")
Treg APOE.CTSZ.TAM
APOE
sp1 = ggplot(data, aes(x = APOE, y = HIF1A))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  xlab("APOE") + 
  ylab('HIF1A') 

sp2 = ggplot(data, aes(x = CTSZ, y = HIF1A))
sp2 = sp2 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  xlab("CTSZ") + 
  ylab('HIF1A') 

sp3 = ggplot(data, aes(x = GLUL, y = HIF1A))
sp3 = sp3 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  xlab("GLUL") + 
  ylab('HIF1A') 


p = plot_grid(sp1, sp2, sp3, nrow = 1)
ggsave("./mmrp.pdf", width = 7.5, height =2.5, device = cairo_pdf, units = 'in')