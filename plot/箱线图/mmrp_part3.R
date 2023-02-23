bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_14", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Pyruvate | Oxaloacetate")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0000621, label='***', alpha=1) + coord_cartesian(ylim = c(0.0000615, 0.00006225))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_28", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
        
            title = "Serine | Methionine")

p=p+rotate_x_text(60)#x轴倾斜角度
p2=p
p2=p2+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0002498, label='***', alpha=1) + coord_cartesian(ylim = c(0.0002474, 0.00025))+
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_40", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Aspartate | Oxaloacetate")

p=p+rotate_x_text(60)#x轴倾斜角度
p3=p
p3=p3+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) + #coord_cartesian(ylim = c(0.00000001, 0.000000016))+
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.00017, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_61", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Putrescine | GABA")

p=p+rotate_x_text(60)#x轴倾斜角度
p4=p
p4=p4+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.04, label='***', alpha=1) +
  theme(legend.position = 'none')

p_all = plot_grid(p1, p2, p3, p4, nrow = 1)
ggsave(filename = "./mmrp_part3.pdf", plot = p_all, device = 'pdf', width = 12, height = 3.5, units = 'in')