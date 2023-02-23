bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"grade"])))]
p=ggboxplot(data, x="grade", y="mean_score", color = "grade",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=grade),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method='wilcox.test', label.x = 0.5)
p1=p1+geom_point(aes(color=grade),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./module1_score_grade.pdf", plot = p1, device = 'pdf', width = 2.5, height = 3.5, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"node_type"])))]
p=ggboxplot(data, x="node_type", y="mean_score", color = "node_type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=node_type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method='wilcox.test', label.x = 0.5)
p1=p1+geom_point(aes(color=node_type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./module1_score_node_type.pdf", plot = p1, device = 'pdf', width = 2.5, height = 3.5, units = 'in')

bioCol=c("#ff3399","#1c79c0","#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"stage_type"])))]
p=ggboxplot(data, x="stage_type", y="mean_score", color = "stage_type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=stage_type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method='wilcox.test', label.x = 0.5)
p1=p1+geom_point(aes(color=stage_type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./module1_score_stage_type.pdf", plot = p1, device = 'pdf', width = 2.5, height = 3.5, units = 'in')

bioCol=c("#ff3399","#1c79c0","#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="type", y="mean_score", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method='wilcox.test', label.x = 0.5)
p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./module2_score_type.pdf", plot = p1, device = 'pdf', width = 2.5, height = 3.5, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"grade_type"])))]
p=ggboxplot(data, x="grade_type", y="mean_score", color = "grade_type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=grade_type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method='wilcox.test', label.x = 0.5)
p1=p1+geom_point(aes(color=grade_type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./module2_score_grade_type.pdf", plot = p1, device = 'pdf', width = 2.5, height = 3.5, units = 'in')
