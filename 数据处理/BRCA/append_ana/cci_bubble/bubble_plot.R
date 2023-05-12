## bubble plot
data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellphonedb/bubble_base.csv", row.names=1)
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
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('gene_expr', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggsave("/public/home/yuwenqi/sc-data/selected/append_ana/cci-bubble/BC/cellphonedb/bubble_bc_cellphonedb.pdf", plot1, width = 10, height = 13, device = cairo_pdf, limitsize=F, units = 'in')
