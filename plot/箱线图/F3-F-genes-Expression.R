bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="M2", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "M2")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.24, label='***', alpha=1) + coord_cartesian(ylim = c(0, 0.25))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="APOE", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
        
            title = "APOE")

p=p+rotate_x_text(60)#x轴倾斜角度
p2=p
p2=p2+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =3.8, label='***', alpha=1) + #coord_cartesian(ylim = c(0.0002, 0.0005))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="CTSZ", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "CTSZ")

p=p+rotate_x_text(60)#x轴倾斜角度
p3=p
p3=p3+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + #coord_cartesian(ylim = c(0.00000001, 0.000000016))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =2.8, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="GPNMB", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "GPNMB")

p=p+rotate_x_text(60)#x轴倾斜角度
p4=p
p4=p4+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0., 3))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =2.4, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="CCL18", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "CCL18")

p=p+rotate_x_text(60)#x轴倾斜角度
p5=p
p5=p5+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0, 2))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =1.6, label='***', alpha=1) + 
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="Anti.inflammatory", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Anti-inflammatory")

p=p+rotate_x_text(60)#x轴倾斜角度
p6=p
p6=p6+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.58, label='***', alpha=1) + #coord_cartesian(ylim = c(0.000180275, 0.000180317))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="CD163", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
        
            title = "CD163")

p=p+rotate_x_text(60)#x轴倾斜角度
p7=p
p7=p7+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y = 2, label='***', alpha=1) + coord_cartesian(ylim = c(0, 2.5))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="CSF1R", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "CSF1R")

p=p+rotate_x_text(60)#x轴倾斜角度
p8=p
p8=p8+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0, 2))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =1.8, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="MSR1", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "MSR1")

p=p+rotate_x_text(60)#x轴倾斜角度
p9=p
p9=p9+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0, 2))+
  annotate("text", x = 'Macro_APOE/CTSZ', y = 1.8, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="cluster", y="PLTP", color = "cluster",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "PLTP")

p=p+rotate_x_text(60)#x轴倾斜角度
p10=p
p10=p10+geom_point(aes(color=cluster),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0, 1.5))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =1.3, label='***', alpha=1) + 
  theme(legend.position = 'none')

p_all = plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2)
ggsave(filename = "./F3-F-genes.pdf", plot = p_all, device = 'pdf', width = 20, height = 7, units = 'in')