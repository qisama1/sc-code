#柱状图
p1 = ggplot(data = data,
       aes(x = Description,y = Log.q.value.,fill = Description),
) +
  geom_bar(stat="identity", width=1, position = position_dodge(width = 1))+
  #geom_text(aes(label=Description), size=2, position = position_stack(vjust = 0))+
  coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size = 8),
        panel.background =element_blank(), axis.title =  element_text(size = 8))+ xlab('') + ylab('Cell Percent')+
  scale_fill_manual(values = c("#E9E681","#F4D03F","#D0AD32","#967E2C","#76D7C4",
                             "#45B39D","#52BE80","#229954","#43856C", "#186A3B",
                             "#F4C9E0","#E6B0AA","#EC7063","#AA6D66","#8F4239",
                             "#C6AED7","#8E74A1","#5F5467","#194684","#39629C",
                             "#607699","#424FD5","#B8B3B1","#737276","#334039") )
# coord_flip() 翻转图表
pdf("/public/home/yuwenqi/sc-data/selected/append_ana/bar_plot/BC/test.pdf", width = 5, height = 2)
p1
dev.off()

axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),