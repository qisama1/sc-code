bioCol=c("#AB82FF", "#20B2AA",'#1c79c0','#ff3399','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"Type"])))]
p=ggboxplot(data, x="Type", y="Macro_APOE.CTSZ", color = "Type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Cell Percent")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=Type),position = position_jitterdodge(),size=1.0) +
    theme(legend.position = 'none') + stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.y = 0.75, label.x = 0.5)  #+ stat_summary(fun='mean', geom="point") #"kruskal.test" "wilcox.test"

ggsave(filename = "./cell_percent.pdf", plot = p1, device = 'pdf', width = 2.6, height = 3, units = 'in')