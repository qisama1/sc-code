p1 = ggplot(data,aes(x=type,y=GLS),fill=type)+
  scale_fill_manual(values=c('#FF1493', '#20B2AA'))+
  geom_boxplot(aes(fill=type),show.legend = F,trim = F)+
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
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 2, label.x = 0.5)  + stat_summary(fun='mean', geom="point")
ggsave("./gls_epi_apoe.pdf", width = 1.5, height =3, device = cairo_pdf, units = 'in')

p2 = ggplot(data,aes(x=type,y=GLS),fill=type)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=type),show.legend = F,trim = F)+
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
  #coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 1.8, label.x = 0.5) + stat_summary(fun='mean', geom="point")
ggsave("./gls_epi_tumor_normal.pdf", width = 1.5, height =3, device = cairo_pdf, units = 'in')

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
  coord_cartesian(ylim = c(0, 2))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 2, label.x = 0.5) + stat_summary(fun='mean', geom="point")
ggsave("./gls_epi_mye.pdf", width = 1.5, height =3, device = cairo_pdf, units = 'in')