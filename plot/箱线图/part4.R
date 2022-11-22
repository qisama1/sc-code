bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_12", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Fumarate | Malate")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0002278, label='***', alpha=1) + #coord_cartesian(ylim = c(0.000180275, 0.000180317))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_15", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
        
            title = "3PD | Serine")

p=p+rotate_x_text(60)#x轴倾斜角度
p2=p
p2=p2+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.00045, label='***', alpha=1) + coord_cartesian(ylim = c(0.0002, 0.0005))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_29", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Cysteine | Pyruvate")

p=p+rotate_x_text(60)#x轴倾斜角度
p3=p
p3=p3+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + #coord_cartesian(ylim = c(0.00000001, 0.000000016))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.00075, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_38", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Aspartate | Asparagine")

p=p+rotate_x_text(60)#x轴倾斜角度
p4=p
p4=p4+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0.0001025, 0.000106))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0001058, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_66", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Citruline+Aspartate | Argininosuccinate")

p=p+rotate_x_text(60)#x轴倾斜角度
p5=p
p5=p5+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + coord_cartesian(ylim = c(0.000215, 0.000225))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0002248, label='***', alpha=1) + 
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_76", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Malate_in")

p=p+rotate_x_text(60)#x轴倾斜角度
p6=p
p6=p6+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.000129, label='***', alpha=1) +
  theme(legend.position = 'none')

p_all = plot_grid(p1, p2, p3, p4, p5, p6, nrow = 1)
ggsave(filename = "./mmrd_part4.pdf", plot = p_all, device = 'pdf', width = 28, height = 3.5, units = 'in')