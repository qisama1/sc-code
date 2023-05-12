data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/ESCA/part3_data.csv")

setwd("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/ESCA/esca_dot_plot3/")
sp1 = ggplot(data, aes(x = CD8Teff_2, y = TAM_C5AR1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("CD8Teff_2") + 
  ylab('TAM_C5AR1') 
ggsave("./CD8Teff_2_TAM_C5AR1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = CD8Teff_2, y = TAM_LGALS9))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("CD8Teff_2") + 
  ylab('TAM_LGALS9') 
ggsave("./CD8Teff_2_TAM_LGALS9.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = CD8Teff_2, y = TAM_CD47))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("CD8Teff_2") + 
  ylab('TAM_CD47') 
ggsave("./CD8Teff_2_TAM_CD47.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

# sp1 = ggplot(data, aes(x = CD8Teff_2, y = TAM_LRP1))
# sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
#   stat_cor(method="spearman", size = 2) + 
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(),
#         #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
#         plot.title = element_text(hjust=0.5))+
#     theme_classic() +
#   theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
#   xlab("CD8Teff_2") + 
#   ylab('TAM_LRP1') 
# ggsave("./CD8Teff_2_TAM_LRP1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = CD8Teff_2, y = TAM_SPP1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("CD8Teff_2") + 
  ylab('TAM_SPP1') 
ggsave("./CD8Teff_2_TAM_SPP1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')
