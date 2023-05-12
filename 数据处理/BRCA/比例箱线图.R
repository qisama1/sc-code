setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_gene_scores/plot/")
data = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/module/cluster_sg_gene_scores/s_celltype_gene_type.csv", row.names=1)
data=melt(data, id.vars=c('type'))

#data$stage <- factor(data$stage,levels = c("Low","High"))

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="variable", y="value", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) 
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave(filename = "./brca_subTME_gene_type.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')

