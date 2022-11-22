## ?????????
sp1 = ggplot(data, aes(x = GLS, y = GLUL))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("GLS") + xlim(-0.02, 4)

sp2 = ggplot(data, aes(x = GLUD1, y = GLUL))
sp2 = sp2 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("GLUD1") + xlim(-0.02, 4)


sp3 = ggplot(data, aes(x = GLUD2, y = GLUL))
sp3 = sp3 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("GLUD2") + xlim(-0.02, 4)

sp4 = ggplot(data, aes(x = OGDH, y = GLUL))
sp4 = sp4 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("OGDH") + xlim(-0.02, 4)

sp5 = ggplot(data, aes(x = SUCLG1, y = GLUL))
sp5 = sp5 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("SUCLG1") + xlim(-0.02, 4)

p4=plot_grid(sp1, sp2,sp3,sp4,sp5,  nrow=1)
ggsave(filename = "./corr/corr_gene.pdf", plot = p4, device = 'pdf', width = 10, height = 2, units = 'in')


