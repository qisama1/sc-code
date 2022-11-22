ecm_myCAF = c("FN1","THBS2","LAMB2","THBS1","COL1A1","CXCL12","IGF1","GRN","JAG1","CXCL14","PLAU","SDC1","SDC2","LRP1")

plot_data = cbind(Fib@meta.data[,'Tissue', drop=F], t(as.matrix(Fib[ecm_myCAF,]@assays$RNA@data)))
colnames(plot_data) = colnames(plot_data) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

## 画箱线图
## 5个一组
p1 = ggplot(plot_data,aes(x=Tissue,y=PLAU),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "PLAU")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)
p2 = ggplot(plot_data,aes(x=Tissue,y=SDC1),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "SDC1")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p3 = ggplot(plot_data,aes(x=Tissue,y=SDC2),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "SDC2")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p4 = ggplot(plot_data,aes(x=Tissue,y=LRP1),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "LRP1")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p5 = ggplot(plot_data,aes(x=Tissue,y=CXCL14),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "CXCL14")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

plot_grid(p1, p2, p3, p4, p5, nrow=1)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/diff_compare/ecm_myCAF/PLAU-SDC1-SDC2-LRP1-CXCL14.pdf", width = 5.2, height =3, device = cairo_pdf, units = 'in')
