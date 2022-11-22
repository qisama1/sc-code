## ?????????
sp1 = ggplot(data, aes(x = M_51, y = GLUL))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_51") + xlim(0, 0.04)

sp2 = ggplot(data, aes(x = M_9, y = GLUL))
sp2 = sp2 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_9") + xlim(0, 0.04)


sp3 = ggplot(data, aes(x = M_54, y = GLUL))
sp3 = sp3 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_54") + xlim(0, 0.04)

sp4 = ggplot(data, aes(x = M_55, y = GLUL))
sp4 = sp4 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_55") + xlim(0, 0.04)

sp5 = ggplot(data, aes(x = M_164, y = GLUL))
sp5 = sp5 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_164") + xlim(0, 0.04)

sp6 = ggplot(data, aes(x = M_166, y = GLUL))
sp6 = sp6 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_166") + xlim(0, 0.04)

sp7 = ggplot(data, aes(x = M_10, y = GLUL))
sp7 = sp7 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("M_10") + xlim(0, 0.04)

p4=plot_grid(sp1, sp2,sp3,sp4,sp5,sp6,sp7,  nrow=1)

ggsave(filename = "./corr/corr_module.pdf", plot = p4, device = 'pdf', width = 14, height = 2, units = 'in')


