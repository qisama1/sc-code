setwd("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/plot/")
data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/module/tpm-compare/tpm_module_stage.csv", row.names=1)
data=melt(data, id.vars=c('stage', 'stage2'))

#data$stage <- factor(data$stage,levels = c("Low","High"))

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"stage"])))]
p=ggboxplot(data, x="variable", y="value", color = "stage",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=stage),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
p1=p1+geom_point(aes(color=stage),position = position_jitterdodge(),size=0.5) 
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave(filename = "./pdac_subTME_stage1&2vs3&4.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"stage2"])))]
p=ggboxplot(data, x="variable", y="value", color = "stage2",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=stage2),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
p1=p1+geom_point(aes(color=stage2),position = position_jitterdodge(),size=0.5) 
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave(filename = "./pdac_subTME_stage1vs2&3&4.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')



