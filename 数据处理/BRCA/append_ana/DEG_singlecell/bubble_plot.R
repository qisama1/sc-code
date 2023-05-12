library(ggplot2)
library(patchwork)
library(dplyr)
## bubble plot
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/DEG_cellphonedb.csv", row.names=1)
data = data[, c('ligand_T', 'ligand_N')]
sel_means=data
sel_pval=data
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = c(colorRampPalette(c("white", '#EE0000'), alpha = TRUE)(n = 400))
plot1 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=mean,color=mean)) +
  scale_color_gradientn('gene_expr', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/DEG_cellphonedb.csv", row.names=1)
data = data[, c('receptor_T', 'receptor_N')]
sel_means=data
sel_pval=data
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = c(colorRampPalette(c("white", '#EE0000'), alpha = TRUE)(n = 400))
plot2 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=mean,color=mean)) +
  scale_color_gradientn('gene_expr', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
remove_y <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)
p <- list(
  plot1,
  plot2 + remove_y + theme(legend.position = "none")
)
plots = wrap_plots(p, nrow = 1) + plot_layout(guides = "collect")
ggsave("/public/home/yuwenqi/sc-data/selected/append_ana/DEG-singlecell/BC/cellphonedb/DEG_cellphonedb.pdf", plots, width = 4, height = 6, device = cairo_pdf, limitsize=F, units = 'in')
