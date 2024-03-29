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
library(ggforce)

bk <- c(seq(0, max(data), (max(data) - 0) / 200))
 
heatmap = pheatmap(data, 
                    #annotation_row = annotation, # ?????????????????????????????????
                    #annotation_col=anno, # ?????????????????????????????????
                    show_colnames = T, # ??????????????????
                    #show_rownames=TRUE,  # ??????????????????
                    fontsize=10, # ????????????
                    breaks = bk,
                    color = c(colorRampPalette(colors=c('#FFFFFF', '#e00000'))(length(bk))), # ?????????????????????
                    #color = colorRampPalette(c("#0c00e6",'#FFFFFF',"#e00000"))(100),
                    annotation_legend=TRUE, # ??????????????????
                    border_color=NA,  # ???????????? NA????????????
                    scale="none",  # ???????????????????????????"row"??????????????????"column"??????????????????"none"?????????
                    cluster_rows = T, # ??????????????????
                    cluster_cols = T, # ??????????????????
                    legend_breaks = c(min(data), 0, max(data)),
                    legend_labels = c(round(min(data), 2),0,round(max(data), 2)),
                    cellheight = 20,
                    cellwidth = 20,
                    angle_col = 45,
                    fontsize_col = 12,
                    fontsize_row =12,
                    #display_numbers = T,
                    #gaps_row = c(3,6,9,12,15)
                    #gaps_row = 1:length(rownames(data)),
                    gaps_col = c(1,2,3)                  
)


bk <- c(seq(min(data), -0.001, (0 - min(data)) / 100),seq(0, max(data), (max(data) - 0) / 100))
library(pheatmap)
heatmap = pheatmap(data, 
                    #annotation_row = annotation, # 
                    #annotation_col=anno, # 
                    show_colnames = T, # 
                    #show_rownames=TRUE,  # 
                    fontsize=10, # 
                    breaks = bk,
                    color = c(colorRampPalette(c("#0c00e6","#FFFFFF"))(length(bk) / 2), colorRampPalette(colors=c('#FFFFFF', '#e00000'))(length(bk) / 2)), # ?????????????????????#color = colorRampPalette(c("#0c00e6",'#FFFFFF',"#e00000"))(100),
                    annotation_legend=TRUE, # 
                    border_color=NA,  #  NA
                    scale="none",  # "row column none
                    cluster_rows = F, # 
                    cluster_cols = F, # 
                    legend_breaks = c(min(data), 0, max(data)),
                    legend_labels = c(round(min(data), 2),0,round(max(data), 2)),
                    cellheight = 20,
                    cellwidth = 20,
                    angle_col = 45,
                    fontsize_col = 12,
                    fontsize_row =12,
                    #display_numbers = T,
                    #gaps_row = c(3,6,9,12,15)
                    #gaps_row = 1:length(rownames(data)),
                    gaps_col = c(1,2,3)
                    
 )