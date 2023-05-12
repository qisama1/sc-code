setwd("/public/home/yuwenqi/sc-data/selected/immune_data/07_SKCM/gene_diff/gene_pair/")
fs = list.files('./', '.csv')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    print(filename)
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('response'))
    boxplot_response(data, paste0("./plot/", filename, ".pdf"))
})


boxplot_response = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"response"])))]
    p=ggboxplot(data, x="variable", y="value", color = "response",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=response),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=response),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}

data = read.csv("/public/home/yuwenqi/sc-data/selected/immune_data/23_NSCLC_GSE126044/23_gene_diff/gene_pair/C3_ITGAX_ITGB2.csv", row.names=1)

p1 = ggplot(data,aes(x=response,y=gene),fill=response)+
  scale_fill_manual(values=c("#ff3399", "#0389ff"))+
  geom_boxplot(aes(fill=response),show.legend = F,trim = F)+ #outlier.shape = NA
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle=60,hjust = 1.1,vjust=1.05),
        plot.title = element_text(hjust=0.5))+
  labs(x=NULL,
       y=NULL,
       title = "response")+
  #theme_minimal() 
  #geom_sina(aes(color=grade), size=0.1,alpha=0.2) 
  #coord_cartesian(ylim = c(0, 2))+
  # wi
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='wilcox.test', label.x = 0.5)  + stat_summary(fun='mean', geom="point")
p1=p1+geom_point(aes(color=response),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave("/public/home/yuwenqi/sc-data/selected/immune_data/23_NSCLC_GSE126044/23_gene_diff/gene_pair/plot/test.pdf", width = 4, height =4, device = cairo_pdf, units = 'in')
