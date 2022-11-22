p1 = ggplot(data,aes(x=cluster,y=GLS),fill=cluster)+
  scale_fill_manual(values=c('#FF1493', '#20B2AA'))+
  geom_boxplot(aes(fill=cluster),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "GLS")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 4))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 3.8, label.x = 0.5) + stat_summary(fun='mean', geom="point")

p2 = ggplot(data,aes(x=cluster,y=GLUD1),fill=cluster)+
  scale_fill_manual(values=c('#FF1493', '#20B2AA'))+
  geom_boxplot(aes(fill=cluster),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "GLUD1")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 1.8, label.x = 0.5) +
  stat_summary(fun='mean', geom="point")



p4 = ggplot(data,aes(x=cluster,y=OGDH),fill=cluster)+
  scale_fill_manual(values=c('#FF1493', '#20B2AA'))+
  geom_boxplot(aes(fill=cluster),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "OGDH")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 1))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 0.8, label.x = 0.5) +
  stat_summary(fun='mean', geom="point")

p5 = ggplot(data,aes(x=cluster,y=SUCLG1),fill=cluster)+
  scale_fill_manual(values=c('#FF1493', '#20B2AA'))+
  geom_boxplot(aes(fill=cluster),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "SUCLG1")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 1.8, label.x = 0.5) +
  stat_summary(fun='mean', geom="point")

plot_grid(p1, p2, p4, p5, nrow=1)
ggsave("./all_genes_epi_mye.pdf", width = 6.5, height =3, device = cairo_pdf, units = 'in')
