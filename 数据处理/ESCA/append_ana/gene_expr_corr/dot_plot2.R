data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/CRC/part2_data.csv")
setwd("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/CRC/crc_dot_plot2/")
sp1 = ggplot(data, aes(x = Macro_APOE, y = Epi_ANXA1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_APOE") + 
  ylab('Epi_ANXA1') 
ggsave("./Macro_APOE_Epi_ANXA1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

# sp1 = ggplot(data, aes(x = Macro_APOE, y = Epi_MDK))
# sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
#   stat_cor(method="spearman", size = 2) + 
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(),
#         #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
#         plot.title = element_text(hjust=0.5))+
#     theme_classic() +
#   theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
#   xlab("Macro_APOE") + 
#   ylab('Epi_MDK') 
# ggsave("./Macro_APOE_Epi_MDK.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = Epi_PLAU))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_APOE") + 
  ylab('Epi_PLAU') 
ggsave("./Macro_APOE_Epi_PLAU.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')
