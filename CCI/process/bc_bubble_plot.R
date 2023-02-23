## bubble plot
sel_means=read.csv("./bc/module1.csv", row.names=1)
sel_pval=read.csv("./bc/module1.csv", row.names=1)
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = colorRampPalette(c("white", "blue", 'yellow', 'red'), alpha = TRUE)(n = 399)

plot1 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

sel_means=read.csv("./bc/module2.csv", row.names=1)
sel_pval=read.csv("./bc/module2.csv", row.names=1)
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = colorRampPalette(c("white", "blue", 'yellow', 'red'), alpha = TRUE)(n = 399)

plot2 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

sel_means=read.csv("./bc/module3.csv", row.names=1)
sel_pval=read.csv("./bc/module3.csv", row.names=1)
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = colorRampPalette(c("white", "blue", 'yellow', 'red'), alpha = TRUE)(n = 399)

plot3 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

sel_means=read.csv("./bc/module4.csv", row.names=1)
sel_pval=read.csv("./bc/module4.csv", row.names=1)
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = colorRampPalette(c("white", "blue", 'yellow', 'red'), alpha = TRUE)(n = 399)

plot4 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

sel_means=read.csv("./bc/others.csv", row.names=1)
sel_pval=read.csv("./bc/others.csv", row.names=1)
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = colorRampPalette(c("white", "blue", 'yellow', 'red'), alpha = TRUE)(n = 399)

plot5 = ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))



ggsave("./bc/cellchat.pdf",plot_grid(plot1, plot2, plot3, plot4, nrow = 1), width = 150, height = 20, device = cairo_pdf, limitsize=F, units = 'in')
ggsave("./bc/cellchat_others.pdf", plot5, width = 250, height = 20, device = cairo_pdf, limitsize=F, units = 'in')
