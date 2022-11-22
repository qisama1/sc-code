#Cycling
Cycling = c("NAMPT","PLXNB2","EFNA1","EFNA3","EFNA4","SDC4","CD44","F11R","DAG1","CXCR4","SDC1","SORT1","JAG1","CD46","NOTCH2","THBS1","GRN","CD47","IGF1R","EGFR","TNFRSF1A","NOTCH3","PLAUR","MDK")

#Cycling
Cycling = c("NAMPT","PLXNB2","EFNA1","EFNA3","EFNA4","SDC4","CD44","F11R","DAG1","CXCR4",
"SDC1","SORT1","JAG1","CD46","NOTCH2","THBS1","GRN","CD47","IGF1R","EGFR","TNFRSF1A","NOTCH3","PLAUR","MDK")

plot_data = cbind(Epi@meta.data[,'Tissue', drop=F], t(as.matrix(Epi[Cycling,]@assays$RNA@data)))
colnames(plot_data) = colnames(plot_data) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

## 画箱线图
## 5个一组
p1 = ggplot(plot_data,aes(x=Tissue,y=TNFRSF1A),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "TNFRSF1A")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)
p2 = ggplot(plot_data,aes(x=Tissue,y=NOTCH3),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "NOTCH3")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p3 = ggplot(plot_data,aes(x=Tissue,y=PLAUR),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "PLAUR")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p4 = ggplot(plot_data,aes(x=Tissue,y=MDK),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "MDK")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

p5 = ggplot(plot_data,aes(x=Tissue,y=EGFR),fill=Tissue)+
  scale_fill_manual(values=c('#20B2AA', '#FF1493'))+
  geom_boxplot(aes(fill=Tissue),show.legend = F)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "EGFR")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  coord_cartesian(ylim = c(0, 5))+
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.y = 5, label.x = 0.5)

plot_grid(p1, p2, p3, p4, p5, nrow=1)
ggsave("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/diff_compare/Cycling/TNFRSF1A-NOTCH3-CD47-PLAUR-MDK.pdf", width = 5.2, height =3, device = cairo_pdf, units = 'in')