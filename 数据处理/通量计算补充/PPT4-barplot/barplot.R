library(ggsci)
library(ggplot2)
library(cowplot)
library(data.table)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(tidyverse)


data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/bc/bc_barplot_Macro_APOE.csv")
## barplot
dat_plot <- data.frame(id = data['X'],
                       t = data['X0'])
# 去掉"HALLMARK_"
library(stringr)
#dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$X0>-0.015, ifelse(dat_plot$X0 >= 0.015 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(X0)
# 变成因子类型
dat_plot$X <- factor(dat_plot$X,levels = dat_plot$X)
# 绘制
library(ggplot2)
#library(ggtheme)
# install.packages("ggprism")

p <- ggplot(data = dat_plot,aes(x = X,y = X0,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-0.015,0.015),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('Matabolism Value') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签
# 此处参考了：https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
# 小于-2的数量
low1 <- dat_plot %>% filter(X0 < -0.015) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter(X0 < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(X0 < 0.015) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = X,y = 0.0001,label = X),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = X,y = 0.0001,label = X),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = X,y = -0.0001,label = X),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = X,y = -0.0001,label = X),
            hjust = 1,color = 'black') # 大于1的为黑色标签
ggsave("/public/home/yuwenqi/sc-data/selected/append_ana/pathway_appends/PPT4-barplot/bc/bc_barplot_Macro_APOE.pdf",p,width = 4,height  = 4)