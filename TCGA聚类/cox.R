library("survival")
#以R包survival中的肺癌数据为例
data("lung")
head(lung)
head(lung)


res<- coxph(Surv(time, status) ~ ph.ecog, data = lung)
summary(res)

# 多个单因素合并
y<- Surv(lung$time, lung$status)
Unicox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  cox<- coxph(FML,data=lung)#
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],2)#指定小数点后数位
  Pvalue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],2)
  CI95 <-round(cox1$conf.int[,4],2)
  CI<- paste0(HR," (",CI5,"-",CI95,")")
  Variable=row.names(data.frame(cox1$coefficients))
  Unicox_model<- data.frame('Variable' = Variable,
                             'HR(95% CI)' = CI,
                             'Pvalue' = Pvalue)
   return(Unicox_model)}  
variable<- c("age","sex",  "ph.karno", "ph.ecog", "wt.loss")#设置变量
Unicox <- lapply(variable, Unicox_model)#循环输入
Unicox<- ldply(Unicox,data.frame)#将list转变为数据框
colnames(Unicox)=c('Variable','HR (95% CI)','Pvalue')
#查看单因素cox表格
View(Unicox)

# 多因素

res <- coxph(Surv(time, status) ~ age + sex + ph.karno + ph.ecog + wt.loss, data =  lung)
mul_cox <- summary(res)
mul_HR<- round(mul_cox$coefficients[,2],2) 
mul_Pvalue<- round(mul_cox$coefficients[,5],4) 
mul_CI5<-round(mul_cox$conf.int[,3],2)
mul_CI95<-round(mul_cox$conf.int[,4],2)
mul_CI<-paste0(mul_HR,' (',mul_CI5,'-',mul_CI95,')')
Variable<-row.names(data.frame(mul_cox$coefficients))
mulcox_res<- data.frame(Variable,mul_HR,mul_CI5,mul_CI95,mul_CI,mul_Pvalue)
colnames(mulcox_res)=c("Variable","HR","CI5","CI95","HR (95% CI)","Pvalue")

# 结果可视化

library(forestplot)
#添加表头
dat=rbind(c("Variable", NA,NA,NA,"HR (95% CI)", "Pvalue"),mulcox_res)
#画图
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
                           "7"= gpar(lty=1,lwd=1.5)), # 在1,2,7行添加横线
           col=fpColors ( box = 'blue ' , #方块颜色 
                          lines = ' black ' ,#置信区间横线颜色
                          zero = "grey" ),#参考线颜色
           lwd.zero=1,#参考线宽度
           lwd.ci=1.5, # 置信区间横线宽度
           lty.ci=7 ,# 置信区间横线类型
           ci.vertices.height=0.1   #置信区间横线两端竖线高度     
)


"/public2022/chenzixi/research/TCGA_data/new-RNAseq/id_name.csv"