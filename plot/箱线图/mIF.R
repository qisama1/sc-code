data=melt(data, id.vars=c("Patients"))

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"variable"])))]     
p=ggboxplot(data, x="variable", y="value", color = "variable",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="          Pearson's correlation 
        coefficient(PCC)",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=variable),position = position_jitterdodge(),size=1.0) +
    theme(legend.position = 'none') + stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.y = 0.9, label.x = 0.5)  #+ stat_summary(fun='mean', geom="point") #"kruskal.test" "wilcox.test"

ggsave(filename = "./pearson.pdf", plot = p1, device = 'pdf', width = 2.8, height = 3.6, units = 'in')


data=melt(data, id.vars=c("Patients"))
bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"variable"])))]     
p=ggboxplot(data, x="variable", y="value", color = "variable",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="Overlap coefficient",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=variable),position = position_jitterdodge(),size=1.0) +
    theme(legend.position = 'none') + stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.y = 0.9, label.x = 0.5)  #+ stat_summary(fun='mean', geom="point") #"kruskal.test" "wilcox.test"

ggsave(filename = "./overlap.pdf", plot = p1, device = 'pdf', width = 2.8, height = 3.6, units = 'in')
