# 散点图
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/PDAC/part1_data.csv")
setwd("/public/home/yuwenqi/sc-data/selected/append_ana/corr-gene-cell/PDAC/pdac_dot_plot/")
CAF_COL1A1, 
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('CAF_COL1A1') 
ggsave("./Ductal_cell_A_1_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('CAF_COL1A1') 
ggsave("./Ductal_cell_A_2_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('CAF_COL1A1') 
ggsave("./Ductal_cell_A_3_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('CAF_COL1A1') 
ggsave("./Ductal_cell_A_4_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('CAF_COL1A1') 
ggsave("./Macro_ACP5_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('CAF_COL1A1') 
ggsave("./Macro_CXCL10_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = CAF_COL1A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('CAF_COL1A1') 
ggsave("./Macro_SPP1_CAF_COL1A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')



CAF_COL1A2, 
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('CAF_COL1A2') 
ggsave("./Ductal_cell_A_1_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('CAF_COL1A2') 
ggsave("./Ductal_cell_A_2_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('CAF_COL1A2') 
ggsave("./Ductal_cell_A_3_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('CAF_COL1A2') 
ggsave("./Ductal_cell_A_4_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('CAF_COL1A2') 
ggsave("./Macro_ACP5_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('CAF_COL1A2') 
ggsave("./Macro_CXCL10_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = CAF_COL1A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('CAF_COL1A2') 
ggsave("./Macro_SPP1_CAF_COL1A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')
CAF_COL6A3, 
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('CAF_COL6A3') 
ggsave("./Ductal_cell_A_1_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('CAF_COL6A3') 
ggsave("./Ductal_cell_A_2_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('CAF_COL6A3') 
ggsave("./Ductal_cell_A_3_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('CAF_COL6A3') 
ggsave("./Ductal_cell_A_4_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('CAF_COL6A3') 
ggsave("./Macro_ACP5_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('CAF_COL6A3') 
ggsave("./Macro_CXCL10_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = CAF_COL6A3))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('CAF_COL6A3') 
ggsave("./Macro_SPP1_CAF_COL6A3.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

CAF_FN1, 
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('CAF_FN1') 
ggsave("./Ductal_cell_A_1_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('CAF_FN1') 
ggsave("./Ductal_cell_A_2_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('CAF_FN1') 
ggsave("./Ductal_cell_A_3_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('CAF_FN1') 
ggsave("./Ductal_cell_A_4_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('CAF_FN1') 
ggsave("./Macro_ACP5_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('CAF_FN1') 
ggsave("./Macro_CXCL10_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = CAF_FN1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('CAF_FN1') 
ggsave("./Macro_SPP1_CAF_FN1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

CAF_THBS2, 
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('CAF_THBS2') 
ggsave("./Ductal_cell_A_1_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('CAF_THBS2') 
ggsave("./Ductal_cell_A_2_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('CAF_THBS2') 
ggsave("./Ductal_cell_A_3_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('CAF_THBS2') 
ggsave("./Ductal_cell_A_4_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('CAF_THBS2') 
ggsave("./Macro_ACP5_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('CAF_THBS2') 
ggsave("./Macro_CXCL10_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = CAF_THBS2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('CAF_THBS2') 
ggsave("./Macro_SPP1_CAF_THBS2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

End_COL4A1, 
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('End_COL4A1') 
ggsave("./Ductal_cell_A_1_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('End_COL4A1') 
ggsave("./Ductal_cell_A_2_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('End_COL4A1') 
ggsave("./Ductal_cell_A_3_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('End_COL4A1') 
ggsave("./Ductal_cell_A_4_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('End_COL4A1') 
ggsave("./Macro_ACP5_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('End_COL4A1') 
ggsave("./Macro_CXCL10_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = End_COL4A1))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('End_COL4A1') 
ggsave("./Macro_SPP1_End_COL4A1.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

End_COL4A2
sp1 = ggplot(data, aes(x = Ductal_cell_A_1, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_1") + 
  ylab('End_COL4A2') 
ggsave("./Ductal_cell_A_1_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_2, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_2") + 
  ylab('End_COL4A2') 
ggsave("./Ductal_cell_A_2_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_3, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_3") + 
  ylab('End_COL4A2') 
ggsave("./Ductal_cell_A_3_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Ductal_cell_A_4, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Ductal_cell_A_4") + 
  ylab('End_COL4A2') 
ggsave("./Ductal_cell_A_4_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_ACP5, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_ACP5") + 
  ylab('End_COL4A2') 
ggsave("./Macro_ACP5_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_CXCL10, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_CXCL10") + 
  ylab('End_COL4A2') 
ggsave("./Macro_CXCL10_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')

sp1 = ggplot(data, aes(x = Macro_SPP1, y = End_COL4A2))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm,  colour = 'black', level = 0.99) +
  stat_cor(method="spearman", size = 2) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
    theme_classic() +
  theme(axis.text = element_text(size = 4), text = element_text(size = 5)) + 
  xlab("Macro_SPP1") + 
  ylab('End_COL4A2') 
ggsave("./Macro_SPP1_End_COL4A2.pdf", width = 1.5, height =1.5, device = cairo_pdf, units = 'in')