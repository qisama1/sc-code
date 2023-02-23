bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="cluster", y="GLS", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))+
  scale_x_discrete(limits=c("Mye", "T/NK", "Epi", "Fib"))
#size width=3,height=6
ggsave(filename = "./GLS_major.pdf", plot = p1, device = 'pdf', width = 6, height = 3, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="cluster", y="GLUD1", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))+
  scale_x_discrete(limits=c("Mye", "T/NK", "Epi", "Fib"))
#size width=3,height=6
ggsave(filename = "./GLUD1_major.pdf", plot = p1, device = 'pdf', width = 6, height = 3, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="cluster", y="GLUD2", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))+
  scale_x_discrete(limits=c("Mye", "T/NK", "Epi", "Fib"))
#size width=3,height=6
ggsave(filename = "./GLUD2_major.pdf", plot = p1, device = 'pdf', width = 6, height = 3, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="cluster", y="OGDH", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))+
  scale_x_discrete(limits=c("Mye", "T/NK", "Epi", "Fib"))
#size width=3,height=6
ggsave(filename = "./OGDH_major.pdf", plot = p1, device = 'pdf', width = 6, height = 3, units = 'in')

bioCol=c("#1c79c0","#ff3399", "#0389ff")
bioCol=bioCol[1:length(levels(factor(data[,"type"])))]
p=ggboxplot(data, x="cluster", y="SUCLG1", color = "type",#??????????????????????????????
            
            xlab="",#x?????????
            
            ylab="",#y?????????
            
            legend.title="",#????????????
            
            palette = bioCol,#??????
            
            width=0.7)

#p=p+rotate_x_text(60)#x???????????????
p1=p+stat_compare_means(aes(group=type),
                        
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        
                        label = "p.signif")
#p1=p1+geom_point(aes(color=type),position = position_jitterdodge(),size=0.5) +
p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))+
  scale_x_discrete(limits=c("Mye", "T/NK", "Epi", "Fib"))
#size width=3,height=6
ggsave(filename = "./SUCLG1_major.pdf", plot = p1, device = 'pdf', width = 6, height = 3, units = 'in')