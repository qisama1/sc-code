library(survivalROC)

#输入数据
#这里我预测五年的生存率
survival_ROC_1<-survivalROC(Stime=data_risk$sur_time, #生存时间，Event time or censoring time for subjects
                          status=data_risk$status, #生存状态,dead or alive
                          marker=data_risk$risk, #风险得分，Predictor or marker value
                          predict.time=1, #预测5年的生存时间
                          method="NNE", #使用KM法进行拟合，默认的方法是method="NNE",
                          span = 0.25 * nrow(ovarian)^(-0.20)
                          )
survival_ROC_AUC_1<-round(survival_ROC_1$AUC,3)#ROC曲线的AUC保留3位小数（文章保留了3位）

survival_ROC_3<-survivalROC(Stime=data_risk$sur_time, #生存时间，Event time or censoring time for subjects
                          status=data_risk$status, #生存状态,dead or alive
                          marker=data_risk$risk, #风险得分，Predictor or marker value
                          predict.time=3, #预测5年的生存时间
                          method="NNE", #使用KM法进行拟合，默认的方法是method="NNE",
                          span = 0.25 * nrow(ovarian)^(-0.20)
                          )
survival_ROC_AUC_3<-round(survival_ROC_3$AUC,3)#ROC曲线的AUC保留3位小数（文章保留了3位）

survival_ROC_5<-survivalROC(Stime=data_risk$sur_time, #生存时间，Event time or censoring time for subjects
                          status=data_risk$status, #生存状态,dead or alive
                          marker=data_risk$risk, #风险得分，Predictor or marker value
                          predict.time=5, #预测5年的生存时间
                          method="NNE", #使用KM法进行拟合，默认的方法是method="NNE",
                          span = 0.25 * nrow(ovarian)^(-0.20)
                          )
survival_ROC_AUC_5<-round(survival_ROC_5$AUC,3)#ROC曲线的AUC保留3位小数（文章保留了3位）

survival_ROC_10<-survivalROC(Stime=data_risk$sur_time, #生存时间，Event time or censoring time for subjects
                          status=data_risk$status, #生存状态,dead or alive
                          marker=data_risk$risk, #风险得分，Predictor or marker value
                          predict.time=10, #预测5年的生存时间
                          method="NNE", #使用KM法进行拟合，默认的方法是method="NNE",
                          span = 0.25 * nrow(ovarian)^(-0.20)
                          )
survival_ROC_AUC_10<-round(survival_ROC_10$AUC,3)
 #画图
#x轴为False positive rate，y轴为True positive rate

time_ROC.res<-data.frame(TP_1year=survival_ROC_1$TP, #获取3年的ROC的TP
                         FP_1year=survival_ROC_1$FP,  #获取3年的ROC的FP
                         TP_3year=survival_ROC_3$TP,  #获取5年的ROC的TP
                         FP_3year=survival_ROC_3$FP, #获取5年的ROC的FP
                         TP_5year=survival_ROC_5$TP, #获取10年的ROC的TP
                         FP_5year=survival_ROC_5$FP,
                         TP_10year=survival_ROC_10$TP, #获取10年的ROC的TP
                         FP_10year=survival_ROC_10$FP) #获取10年的ROC的FP


TimeROC_plot<-ggplot()+
  geom_line(data=time_ROC.res,aes(x=FP_1year,y=TP_1year),size=1,color="red")+
  geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=1,color="blue")+
  geom_line(data=time_ROC.res,aes(x=FP_5year,y=TP_5year),size=1,color="black")+
  geom_line(data=time_ROC.res,aes(x=FP_10year,y=TP_10year),size=1,color="green")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2 #添加虚线
            )+
  theme_bw()+
  labs(x="False positive rate",y="True positive rate")+
  annotate("text",x=0.75,y=0.3,size=4.5,
           label=paste0("AUC at 1 years = ",survival_ROC_AUC_1),color="red") +
  annotate("text",x=0.75,y=0.25,size=4.5,
           label=paste0("AUC at 3 years = ",survival_ROC_AUC_3),color="blue") +
  annotate("text",x=0.75,y=0.20,size=4.5,
           label=paste0("AUC at 5 years = ",survival_ROC_AUC_5),color="black")+
  annotate("text",x=0.75,y=0.15,size=4.5,
           label=paste0("AUC at 10 years = ",survival_ROC_AUC_10),color="green")+
  theme(axis.text=element_text(face="bold", size=11,  color="black"),#加粗刻度标签
        axis.title=element_text(face="bold", size=14, color="black"))  #加粗xy轴标签名字

ggsave(filename = "./roc.pdf", plot = TimeROC_plot, device = 'pdf', width = 15, height = 15, units = 'cm')