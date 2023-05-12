#柱状图
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/bar_plot/BC/subTME1.csv")
p1 = ggplot(data = data,
       aes(x = Description,y = Log.q.value.,fill = Description),
) +
  geom_bar(stat="identity", width=1, position = position_dodge(width = 1))+
  #geom_text(aes(label=Description), size=2, position = position_stack(vjust = 0))+
  coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size = 8),
        panel.background =element_blank(), axis.title =  element_text(size = 8))+ xlab('') + ylab('subTME1')+
  scale_fill_manual(values = c("#E9E681","#F4D03F","#39629C","#967E2C","#76D7C4",
                             "#45B39D","#EC7063","#229954","#424FD5", "#B8B3B1",
                             "#F4C9E0","#E6B0AA","#EC7063","#AA6D66","#8F4239",
                             "#C6AED7","#8E74A1","#5F5467","#194684","#39629C",
                             "#607699","#424FD5","#B8B3B1","#737276","#334039") )

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/bar_plot/BC/subTME2.csv")
p2 = ggplot(data = data,
       aes(x = Description,y = Log.q.value.,fill = Description),
) +
  geom_bar(stat="identity", width=1, position = position_dodge(width = 1))+
  #geom_text(aes(label=Description), size=2, position = position_stack(vjust = 0))+
  coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size = 8),
        panel.background =element_blank(), axis.title =  element_text(size = 8))+ xlab('') + ylab('subTME2')+
  scale_fill_manual(values = c("#E9E681","#F4D03F","#39629C","#967E2C","#76D7C4",
                             "#45B39D","#EC7063","#229954","#424FD5", "#B8B3B1",
                             "#F4C9E0","#E6B0AA","#EC7063","#AA6D66","#8F4239",
                             "#C6AED7","#8E74A1","#5F5467","#194684","#39629C",
                             "#607699","#424FD5","#B8B3B1","#737276","#334039") )

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/bar_plot/BC/subTME4.csv")
p3 = ggplot(data = data,
       aes(x = Description,y = Log.q.value.,fill = Description),
) +
  geom_bar(stat="identity", width=1, position = position_dodge(width = 1))+
  #geom_text(aes(label=Description), size=2, position = position_stack(vjust = 0))+
  coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size = 8),
        panel.background =element_blank(), axis.title =  element_text(size = 8))+ xlab('') + ylab('subTME4')+
  scale_fill_manual(values = c("#E9E681","#F4D03F","#39629C","#967E2C","#76D7C4",
                             "#45B39D","#EC7063","#229954","#424FD5", "#B8B3B1",
                             "#F4C9E0","#E6B0AA","#EC7063","#AA6D66","#8F4239",
                             "#C6AED7","#8E74A1","#5F5467","#194684","#39629C",
                             "#607699","#424FD5","#B8B3B1","#737276","#334039") )

p <- list(
  p1, p2, p3
)
plots = wrap_plots(p, nrow = 3) + plot_layout(guides = "collect")
ggsave("/public/home/yuwenqi/sc-data/selected/append_ana/bar_plot/BC/bc_subTME_barplot.pdf", plots, width = 4, height = 6, device = cairo_pdf, limitsize=F, units = 'in')
