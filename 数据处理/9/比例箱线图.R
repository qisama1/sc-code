data = read.csv("./major.csv", row.names=1)
data=melt(data, id.vars=c('Grade', 'Stage', 'stage_node'))

#data$stage <- factor(data$stage,levels = c("Low","High"))

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"Grade"])))]
p=ggboxplot(data, x="variable", y="value", color = "Grade",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=Grade),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave(filename = "./major_Grade.pdf", plot = p1, device = 'pdf', width = 24, height = 3.5, units = 'in')



#data = read.csv("./major_tumor.csv", row.names=1)
#data=melt(data, id.vars=c('type', 'stage'))

#data$stage <- factor(data$stage,levels = c("Low","High"))

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"Stage"])))]
p=ggboxplot(data, x="variable", y="value", color = "Stage",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=Stage),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave(filename = "./major_Stage.pdf", plot = p1, device = 'pdf', width = 24, height = 3.5, units = 'in')



#data = read.csv("./major_tumor.csv", row.names=1)
#data=melt(data, id.vars=c('type', 'stage', 'grade'))

#data$grade <- factor(data$grade,levels = c("Low","High"))

#data$stage <- factor(data$stage,levels = c("Low","High"))

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"stage_node"])))]
p=ggboxplot(data, x="variable", y="value", color = "stage_node",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=stage_node),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

ggsave(filename = "./major_stage_node.pdf", plot = p1, device = 'pdf', width = 24, height = 3.5, units = 'in')

