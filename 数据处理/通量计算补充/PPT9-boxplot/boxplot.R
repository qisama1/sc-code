


p1 = ggplot(data,aes(x=type,y=SLC1A3),fill=type)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_violin(aes(fill=type),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "SLC1A3")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 1))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 2.5, label.x = 0.5)
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT8-boxplot/bc/SLC1A3_type.pdf", plot = p1, device = 'pdf', width = 1.5, height = 2, units = 'in')

p2 = ggplot(data,aes(x=type,y=GLUL),fill=type)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_violin(aes(fill=type),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "GLUL")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 1))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 3.7, label.x = 0.5)

p3 = ggplot(data,aes(x=type,y=M_91),fill=type)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_violin(aes(fill=type),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "M91")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 1))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 3.5, label.x = 0.5)

p4 = ggplot(data,aes(x=type,y=M_48),fill=type)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=type),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "M48")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 0.1))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 0.25, label.x = 0.5)


p5 = ggplot(data,aes(x=type,y=Glutamine),fill=type)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_violin(aes(fill=type),show.legend = F,trim = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "Glutamine")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 0.3))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 0.25, label.x = 0.5)

p6 = plot_grid(p1, p2, p3, p4, p5, nrow = 1)
ggsave(filename = "mmrp.pdf", plot = p6, device = 'pdf', width = 7, height = 2, units = 'in')

