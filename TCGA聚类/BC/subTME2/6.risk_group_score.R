data_risk = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/risk2/all_risk.csv", row.names=1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/os/os_res.csv")

os_res$days_to_last_followup[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
 
os_res=os_res[os_res$vital_status %in% c('Alive','Dead'),]
os_res$days_to_last_followup[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
os_res$os = ifelse(os_res$vital_status =='Dead', 1, 0)
os_res = os_res[os_res$days > 0,]
#tcga_data = tcga_data[os_res$Sample.ID, ]

data_risk$sur_time = os_res$days / 365
data_risk$status <- os_res$os
data_risk$id <- row.names(data_risk)
k = length(data_risk)
data_risk=arrange(data_risk, risk)
row.names(data_risk) <- data_risk$id
data_risk <- data_risk[,1:k]
q <- length(row.names(data_risk))
s <- 100/q
i <- 1
precent <- 1
while (i<=q){
  precent[i] <- i*s;
  i <- i+1
}



data_risk$precent <- precent
dat1 <- data_risk[data_risk$precent <= 50,]
dat2 <- data_risk[data_risk$precent > 50,]
dat3 <- data_risk[data_risk$status==0,]
dat4 <- data_risk[data_risk$status==1,]
data_risk <- within(data_risk,{   ##gourp precent data
  risk_group <- NA
  risk_group[precent <= 50] <- 1
  risk_group[precent >50] <- 0
})

setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/risk_group_score/")
tiff(file="KM_risk_score.tiff",width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=300)
fit <- survfit(Surv(sur_time, status) ~ risk_group, data = data_risk)
diff <- survdiff(Surv(sur_time,status) ~ risk_group,data = data_risk,rho = 0)

res <- ggsurvplot(fit,data = data_risk,
                    legend.title = "risk score",
                    legend.labs = c("high risk score","low risk score"),
                    pval = TRUE, ##是否显示P值，TRUE可以换为自己想要表达的字符串
                    pval.method = TRUE, # 是否展示P值的检验方法
                    conf.int = TRUE, ##是否展示曲线的置信区间
                    cont.int.style = 'step', ## 自定义置信区间的style
                    xlab = 'Time in years', ## 自定义X轴的label
                    ggtheme = theme_classic(), # 设置主题(自定义ggplot2中的主题)
                    risk.table = T, ## 所有可用的value是TRUE, FALSE, 'absolute', 'percentage', 'nrisk_cumcensor', 'nrisk_cumevents' 
                    risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
                    risk.table.y.text = F, ## 是否展示分层变量的文本，如果不展示，则用颜色条进行展示
                    break.time.by = 2,
                    #surv.median.line = "hv", # 指出中位生存率1F78B4
                    palette = c("#E31A1C","#1F78B4") ## 自定义曲线的颜色
                    ) # Reverse plot
res$plot <- res$plot + labs(title = "Survival Curve of Risk score",x="Time (years)")+
  theme(plot.title = element_text(hjust = 0.5))
res
dev.off()
os_res$os = ifelse(os_res$vital_status == 1, "Dead", "Alive")
data_risk$os = os_res$os

data_risk$status <- as.factor(data_risk$status)
tiff(file="plot_risk.tiff",width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=300)
layout_1 <- grid.layout(nrow = 3, ncol = 1, heights = c(1,4.5,4.5)) 
pushViewport(viewport(layout = layout_1))  
grid.text("Risk score(AP) for muti-cox analysis", x = 0.5, y = 0.95, gp = gpar(col = "black", fontsize = 15))  # ???ӻ???????

plot1 <- ggplot(data_risk,aes(x=precent,y=sur_time,color=os)) +
  geom_point(cex=0.8)+
  geom_point(data=dat4,aes(x=precent,y=sur_time),colour="#E31A1C",cex=1) +
  geom_point(data=dat3,aes(x=precent,y=sur_time),colour="#1F78B4",cex=1) +
  labs(y = "Survival time(years)",x="Risk score percent")+
  theme(axis.title.y = element_text(size=12),legend.position = c(0.925,0.8))+
  scale_color_manual(values = c("Dead" = "#E31A1C", "Alive" = "#1F78B4"))+
  geom_vline(aes(xintercept=50),linetype="dashed")

plot2 <- ggplot(data_risk,aes(x=precent,y=risk)) +
  geom_point(cex=0.5)+
  geom_point(data=dat2,aes(x=precent,y=risk),colour="#E31A1C",cex=1) +
  geom_point(data=dat1,aes(x=precent,y=risk),colour="#1F78B4",cex=1) +
  labs(y = "Risk score",x="Risk score percent")+
  theme(axis.title.y = element_text(size=12))+
  geom_vline(aes(xintercept=50),linetype="dashed")+
  geom_hline(aes(yintercept=median(data_risk$risk)),linetype="dashed")

print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(plot1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
