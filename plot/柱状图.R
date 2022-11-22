#柱状图
ggplot(data = data,
       aes(x = DataSet,y = Percent,fill = DataSet),
) +
  geom_bar(stat="identity", width=0.8, position = position_dodge(0.2))+
  coord_flip() +
  theme(legend.position = "none", axis.text = element_text(color='black'),
        panel.background =element_blank())+ xlab('') + ylab('Cell Percent') +
  scale_color_manual(values = c('#9A32CD','#5CACEE','#CD6090', '#D2691E'))+
  scale_fill_manual(values = c('#9A32CD','#5CACEE','#CD6090', '#D2691E'))
# coord_flip() 翻转图表