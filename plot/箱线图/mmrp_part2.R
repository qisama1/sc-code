bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_20", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Glycine | Glyoxylate")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.00017190, label='***', alpha=1) + coord_cartesian(ylim = c(0.000171890, 0.000171893))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_43", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
        
            title = "Histidine | B-Alanine")

p=p+rotate_x_text(60)#x轴倾斜角度
p2=p
p2=p2+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0000228, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_51", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Glutamate | 2OG")

p=p+rotate_x_text(60)#x轴倾斜角度
p3=p
p3=p3+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0.00000002, 0.00000005))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.000000049, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_60", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "lysine | Acetyl-CoA")

p=p+rotate_x_text(60)#x轴倾斜角度
p4=p
p4=p4+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.000057, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_74", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "2OG_in | 2OG")

p=p+rotate_x_text(60)#x轴倾斜角度
p5=p
p5=p5+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + #coord_cartesian(ylim = c(0.000029, 0.000035))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.000044, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_106", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Glucose | Glucose-6-phosphate")

p=p+rotate_x_text(60)#x轴倾斜角度
p6=p
p6=p6+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.00015, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_153", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,

            title = "UMP | CDP")

p=p+rotate_x_text(60)#x轴倾斜角度
p7=p
p7=p7+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0.00003538, 0.0000357))+
  #annotate("text", x = 'Macro_APOE/CTSZ', y =0.00000148, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_155", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "UTP | CDP")

p=p+rotate_x_text(60)#x轴倾斜角度
p8=p
p8=p8+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0.000025, 0.0000251))+
  #annotate("text", x = 'Macro_APOE/CTSZ', y =0.0000739, label='***', alpha=1) + coord_cartesian(ylim = c(0.000072, 0.000074))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_156", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "CDP | Cytidine")

p=p+rotate_x_text(60)#x轴倾斜角度
p9=p
p9=p9+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0.0001752, 0.0001758))+
  #annotate("text", x = 'Macro_APOE/CTSZ', y =0.0000739, label='***', alpha=1) + coord_cartesian(ylim = c(0.000072, 0.000074))+
  theme(legend.position = 'none')


p_all = plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 1)
ggsave(filename = "./mmrp_part2.pdf", plot = p_all, device = 'pdf', width = 28, height = 3.5, units = 'in')