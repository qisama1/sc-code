sp1 = ggplot(data, aes(x = Treg, y = Percent_Mye))
sp1 = sp1 + geom_point(colour='#F08080') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) 
  ylab("xxx")


sp2 = ggplot(data, aes(x = M2, y = M_48))
sp2 = sp2 + geom_point(colour='#B0E2FF') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylim(0, 0.2)

sp3 = ggplot(data, aes(x = Pro.inflammatory, y = M_48))
sp3 = sp3 + geom_point(colour='#B0E2FF') + stat_smooth(method = lm, se = F, colour = 'black') +
  stat_cor(method="spearman") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5)) +
  ylim(0, 0.2)
p4=plot_grid(sp1, sp2,sp3,  nrow=1)

ggsave(filename = "./ppt19-1/mmrp.pdf", plot = p4, device = 'pdf', width = 6.5, height = 2, units = 'in')