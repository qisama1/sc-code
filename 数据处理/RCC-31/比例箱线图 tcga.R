setwd("/public/home/yuwenqi/sc-data/selected/immune_data/32_RCC/32_RCC_subTME_diff/")
fs = list.files('./', '_response.csv')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    print(filename)
    data = read.csv(x, row.names = 1)
    data=melt(data, id.vars=c('response'))
    boxplot_response(data, paste0("./plot/", filename, ".pdf"))
})

boxplot_response = function(data, filename) {
    bioCol=c("#1c79c0","#ff3399", "#0389ff")
    bioCol=bioCol[1:length(levels(factor(data[,"response"])))]
    p=ggboxplot(data, x="variable", y="value", color = "response",#??????????????????????????????
                
                xlab="",#x?????????
                
                ylab="",#y?????????
                
                legend.title="",#????????????
                
                palette = bioCol,#??????
                
                width=0.6)

    p=p+rotate_x_text(60)#x???????????????
    p1=p+stat_compare_means(aes(group=response),
                            
                            symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method='wilcox.test', label = "p.signif")
    p1=p1+geom_point(aes(color=response),position = position_jitterdodge(),size=0.5) 
    p1=p1+theme(legend.position = c(0.8, 0.9), legend.background = element_blank(), legend.key = element_blank(),
                panel.border = element_rect(linetype = 'solid', size = 0.3,fill = NA))

    ggsave(filename = filename, plot = p1, device = 'pdf', width = 3.5, height = 3.5, units = 'in')

}
