## ?????????
sp1 = ggplot(data, aes(x = X2OG, y = GLUL))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("2OG") + xlim(-0.02, 0.2)

sp2 = ggplot(data, aes(x = Succinyl.CoA, y = GLUL))
sp2 = sp2 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylab("GLUL") + xlab("Succinyl-CoA") + xlim(-0.02, 0.2)

p4=plot_grid(sp1, sp2,  nrow=1)

ggsave(filename = "./corr/corr_metabolism.pdf", plot = p4, device = 'pdf', width = 4, height = 2, units = 'in')


