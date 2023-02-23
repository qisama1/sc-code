## bc
data=melt(data, id.vars=c("orig.ident", 'type', 'module1'))

data$grade <- factor(data$grade,levels = c("low","high"))

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"grade"])))]
p=ggboxplot(data, x="variable", y="value", color = "grade",#??????????????????????????????

            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=grade),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method = 'wilcox.test')
p1=p1+geom_point(aes(color=grade),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
#size width=3,height=6
ggsave(filename = "./score_tumor_normal.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')


bioCol=c("#ff3399","#1c79c0", "#0389ff")
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
                        
                        label = "p.signif", method = 'wilcox.test')
p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./score_grade_1&2_3.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')

bioCol=c("#ff3399","#1c79c0", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type2"])))]
p=ggboxplot(data, x="variable", y="value", color = "type2",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type2),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method = 'kruskal.test')
p1=p1+geom_point(aes(color=type2),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./score_grade_1_2_3.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')

## crc
data=melt(data, id.vars=c("pid", 'type', 'grade', 'node_type', 'stage_type', 'module1'))
bioCol=c("#ff3399","#1c79c0", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"stage_type"])))]
p=ggboxplot(data, x="variable", y="value", color = "stage_type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=stage_type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif" , method = 'wilcox.test')
p1=p1+geom_point(aes(color=stage_type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./score_stage_type.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')


bioCol=c("#ff3399","#1c79c0", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"grade"])))]
p=ggboxplot(data, x="variable", y="value", color = "grade",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=grade),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method = 'wilcox.test')
p1=p1+geom_point(aes(color=grade),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./score_grade.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')



bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"node_type"])))]
p=ggboxplot(data, x="variable", y="value", color = "node_type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=node_type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif", method = 'wilcox.test')
p1=p1+geom_point(aes(color=node_type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./score_node_type.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')


data=melt(data, id.vars=c("pid", 'type', 'module1'))
bioCol=c("#ff3399","#1c79c0", "#0389ff")
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
                        
                        label = "p.signif", method = 'wilcox.test')
p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5)
#p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            #panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))
ggsave(filename = "./crc_score_tumor_normal.pdf", plot = p1, device = 'pdf', width = 6, height = 3.5, units = 'in')
