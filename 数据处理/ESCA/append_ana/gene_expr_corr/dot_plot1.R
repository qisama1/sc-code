# 散点图
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/CRC/part1_data.csv")
setwd("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/CRC/crc_dot_plot/")
CAF_COL1A1, CAF_COL1A2, CAF_COL6A3, CAF_FN1, CAF_THBS2, End_COL4A1, End_COL4A2

sp1 = ggplot(data, aes(x = Cancer_Malig, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Cancer_Malig") + 
  ylab('CAF_COL1A1') 
ggsave("./Cancer_Malig_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = CAF_COL1A1))
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
  ylab('CAF_COL1A1') 
ggsave("./Macro_APOE_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Cancer_Malig, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Cancer_Malig") + 
  ylab('CAF_COL1A2') 
ggsave("./Cancer_Malig_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = CAF_COL1A2))
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
  ylab('CAF_COL1A2') 
ggsave("./Macro_APOE_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Cancer_Malig, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Cancer_Malig") + 
  ylab('CAF_COL6A3') 
ggsave("./Cancer_Malig_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = CAF_COL6A3))
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
  ylab('CAF_COL6A3') 
ggsave("./Macro_APOE_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

# sp1 = ggplot(data, aes(x = Cancer_Malig, y = CAF_FN1))
# sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
#   stat_cor(method="spearman", size = 2) + 
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(),
#         #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
#         plot.title = element_text(hjust=0.5))+
#     theme_classic() +
#   theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
#   xlab("Cancer_Malig") + 
#   ylab('CAF_FN1') 
# ggsave("./Cancer_Malig_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = CAF_FN1))
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
  ylab('CAF_FN1') 
ggsave("./Macro_APOE_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

# sp1 = ggplot(data, aes(x = Cancer_Malig, y = CAF_THBS2))
# sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
#   stat_cor(method="spearman", size = 2) + 
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(),
#         #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
#         plot.title = element_text(hjust=0.5))+
#     theme_classic() +
#   theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
#   xlab("Cancer_Malig") + 
#   ylab('CAF_THBS2') 
# ggsave("./Cancer_Malig_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

# sp1 = ggplot(data, aes(x = Macro_APOE, y = CAF_THBS2))
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
#   ylab('CAF_THBS2') 
# ggsave("./Macro_APOE_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Cancer_Malig, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Cancer_Malig") + 
  ylab('End_COL4A1') 
ggsave("./Cancer_Malig_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = End_COL4A1))
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
  ylab('End_COL4A1') 
ggsave("./Macro_APOE_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Cancer_Malig, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Cancer_Malig") + 
  ylab('End_COL4A2') 
ggsave("./Cancer_Malig_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_APOE, y = End_COL4A2))
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
  ylab('End_COL4A2') 
ggsave("./Macro_APOE_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')