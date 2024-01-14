data_risk = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/risk2/all_risk.csv", row.names=1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/COAD/os/os_meta.csv")
os_res$days_to_last_followup[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
 
os_res=os_res[os_res$vital_status %in% c('Alive','Dead'),]
os_res$days_to_last_followup[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
os_res$os = ifelse(os_res$vital_status =='Dead', 1, 0)
os_res$sex = ifelse(os_res$gender =='FEMALE', 1, 0)
os_res$age = os_res$age_at_initial_pathologic_diagnosis
os_res = os_res[os_res$days > 0,]

data = cbind(data_risk, os_res)
res <- coxph(Surv(days, os) ~ age + sex  + risk, data =  data)
mul_cox <- summary(res)
mul_HR<- round(mul_cox$coefficients[,2],2) 
mul_Pvalue<- round(mul_cox$coefficients[,5],4) 
mul_CI5<-round(mul_cox$conf.int[,3],2)
mul_CI95<-round(mul_cox$conf.int[,4],2)
mul_CI<-paste0(mul_HR,' (',mul_CI5,'-',mul_CI95,')')
Variable<-row.names(data.frame(mul_cox$coefficients))
mulcox_res<- data.frame(Variable,mul_HR,mul_CI5,mul_CI95,mul_CI,mul_Pvalue)
colnames(mulcox_res)=c("Variable","HR","CI5","CI95","HR (95% CI)","Pvalue")

library(forestplot)
#添加表头
dat=rbind(c("Variable", NA,NA,NA,"HR (95% CI)", "Pvalue"),mulcox_res)
#画图
setwd("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/forestplot/")
pdf("./riskscore_age_sex.pdf",onefile=FALSE, width = 5, height = 2)
forestplot(dat[,c(1,5,6)], #显示表格的第1，5，6列内容
           mean=dat[,2],   #第2列为HR，变成森林图的方块
           lower=dat[,3], upper=dat[,4], #第3列为5%CI，第4列为95%CI，将化作线段，穿过方块
           zero=1,            #零线或参考线位置为HR=1
           boxsize=0.2,       #设置方块大小
           graph.pos=3,#将森林图插在第3列
           xticks=c(0,1,2,3) ,# 设置横轴数字
           txt_gp=fpTxtGp (
             label=gpar(cex=0.8) ,ticks=gpar(cex=0.6)
            ),#调整字体
           hrzl_lines=list("1" = gpar(lty=1,lwd=1.5),
                           "2" = gpar(lty=1,lwd=1.5),
                           "5"= gpar(lty=1,lwd=1.5)), # 在1,2,7行添加横线
           col=fpColors ( box = 'blue ' , #方块颜色 
                          lines = ' black ' ,#置信区间横线颜色
                          zero = "grey" ),#参考线颜色
           lwd.zero=1,#参考线宽度
           lwd.ci=1.5, # 置信区间横线宽度
           lty.ci=7 ,# 置信区间横线类型
           ci.vertices.height=0.1   #置信区间横线两端竖线高度     
)
dev.off()