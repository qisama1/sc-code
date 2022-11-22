bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_16", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Serine | Pyruvate")

p=p+rotate_x_text(60)#x轴倾斜角度
p1=p
p1=p1+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.004, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_27", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
        
            title = "Serine | 2-Aminoacrylate")

p=p+rotate_x_text(60)#x轴倾斜角度
p2=p
p2=p2+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0000004, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_32", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "2-Aminoacrylate | Propanoyl-CoA")

p=p+rotate_x_text(60)#x轴倾斜角度
p3=p
p3=p3+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.02, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_48", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            title = "Glutamate | Glutamine")

p=p+rotate_x_text(60)#x轴倾斜角度
p4=p
p4=p4+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.03, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_85", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Aspartate_in")

p=p+rotate_x_text(60)#x轴倾斜角度
p5=p
p5=p5+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.004, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_91", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "glutamate_in")

p=p+rotate_x_text(60)#x轴倾斜角度
p6=p
p6=p6+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.008, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_103", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,

            title = "Proline_in")

p=p+rotate_x_text(60)#x轴倾斜角度
p7=p
p7=p7+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.00000148, label='***', alpha=1) +
  theme(legend.position = 'none')

bioCol=c("#1c79c0", "#ff3399",'#FFDEAD','#7CFC00','#20B2AA', '#006400', '#AB82FF')
bioCol=bioCol[1:length(levels(factor(data[,"cell_type"])))]
p=ggboxplot(data, x="cell_type", y="M_105", color = "cell_type",#根据分型结果定义颜色
            
            xlab="",#x轴名称
            
            ylab="",#y轴名称
            
            legend.title="",#图例标题 
            
            palette = bioCol,#颜色
            
            width=0.4,
            
            title = "Fatty Acid_in")

p=p+rotate_x_text(60)#x轴倾斜角度
p8=p
p8=p8+geom_point(aes(color=cell_type),position = position_jitterdodge(),size=1.0) +
  scale_x_discrete(limits=c("Macro_APOE/CTSZ", "Macro_CXCL8", "Macro_CCL19", "Macro_FTL", "Macro_Pro","Mono")) +
  annotate("text", x = 'Macro_APOE/CTSZ', y =0.0000739, label='***', alpha=1) + coord_cartesian(ylim = c(0.000072, 0.000074))+
  theme(legend.position = 'none')

p_all = plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 1)
ggsave(filename = "./mmrd_part1.pdf", plot = p_all, device = 'pdf', width = 25, height = 3.5, units = 'in')