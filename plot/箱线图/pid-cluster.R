bioCol=c("#AB82FF", "#ff3399",'#1c79c0','#20B2AA','#20B2AA', '#006400', '#AB82FF', '#8B0000')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="CXCL16", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CXCL16")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =2.9, label='***', alpha=1) + #coord_cartesian(ylim = c(0.000072, 0.000074))+
  theme(legend.position = 'none')
ggsave(filename = "./CXCL16.pdf", plot = p1, device = 'pdf', width = 3.5, height = 4, units = 'in')

bioCol=c("#AB82FF", "#ff3399",'#1c79c0','#20B2AA', '#006400', '#8B0000', '#708090', '#4F94CD', '#473C8B')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="CXCR6", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "CXCL16")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  stat_compare_means(symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif", method='kruskal.test', label.y = 0.95, label.x = 0.5) +
  #scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  #annotate("text", x = 'Macro_APOE/CTSZ', y =2.9, label='***', alpha=1) + #coord_cartesian(ylim = c(0.000072, 0.000074))+
  theme(legend.position = 'none')
ggsave(filename = "./CXCR6.pdf", plot = p1, device = 'pdf', width = 3.5, height = 4, units = 'in')