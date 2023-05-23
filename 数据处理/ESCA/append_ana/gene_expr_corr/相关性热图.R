data = read.csv("/public/home/yuwenqi/sc-data/selected/append_ana/5-15-corr-plot/crc/crc_mye.csv", row.names=1)
data = as.matrix(data)
#data = data[, c('FPR3', 'PLAUR' , 'LRP1','CCL3', 'CCL4' , 'C5AR1')]

pdf("/public/home/yuwenqi/sc-data/selected/append_ana/5-15-corr-plot/crc/crc_mye_corr.pdf", width=5.2, height=4.8)
col = c(colorRampPalette(c("#0c00e6","#FFFFFF"))(100), colorRampPalette(colors=c('#EEA9B8', '#ff0c00'))(100))
p1=corrplot(data, method = 'color', order = "hclust", tl.col = "black", tl.srt = 45, col = col, tl.cex = 0.8)
dev.off()