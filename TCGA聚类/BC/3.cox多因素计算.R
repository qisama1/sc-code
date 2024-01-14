

library("survival")
#以R包survival中的肺癌数据为例
data = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/coad_tcga.csv", row.names = 1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/COAD/os/os_meta.csv")

# res <- coxph(Surv(time, status) ~ ., data =  lung)
# mul_cox <- summary(res)
# mul_HR<- round(mul_cox$coefficients[,2],2) 
# mul_Pvalue<- round(mul_cox$coefficients[,5],4) 
# mul_CI5<-round(mul_cox$conf.int[,3],2)
# mul_CI95<-round(mul_cox$conf.int[,4],2)
# mul_CI<-paste0(mul_HR,' (',mul_CI5,'-',mul_CI95,')')
# Variable<-row.names(data.frame(mul_cox$coefficients))
# mulcox_res<- data.frame(Variable,mul_HR,mul_CI5,mul_CI95,mul_CI,mul_Pvalue)
# colnames(mulcox_res)=c("Variable","HR","CI5","CI95","HR (95% CI)","Pvalue")

# # 结果可视化

# library(forestplot)
# #添加表头
# dat=rbind(c("Variable", NA,NA,NA,"HR (95% CI)", "Pvalue"),mulcox_res)
# #画图
# forestplot(dat[,c(1,5,6)], #显示表格的第1，5，6列内容
#            mean=dat[,2],   #第2列为HR，变成森林图的方块
#            lower=dat[,3], upper=dat[,4], #第3列为5%CI，第4列为95%CI，将化作线段，穿过方块
#            zero=1,            #零线或参考线位置为HR=1
#            boxsize=0.2,       #设置方块大小
#            graph.pos=3,#将森林图插在第3列
#            xticks=c(0,1,2,3) ,# 设置横轴数字
#            txt_gp=fpTxtGp (
#              label=gpar(cex=0.8) ,ticks=gpar(cex=0.6)
#             ),#调整字体
#            hrzl_lines=list("1" = gpar(lty=1,lwd=1.5),
#                            "2" = gpar(lty=1,lwd=1.5),
#                            "7"= gpar(lty=1,lwd=1.5)), # 在1,2,7行添加横线
#            col=fpColors ( box = 'blue ' , #方块颜色 
#                           lines = ' black ' ,#置信区间横线颜色
#                           zero = "grey" ),#参考线颜色
#            lwd.zero=1,#参考线宽度
#            lwd.ci=1.5, # 置信区间横线宽度
#            lty.ci=7 ,# 置信区间横线类型
#            ci.vertices.height=0.1   #置信区间横线两端竖线高度     
# )

cox_cal = function(os_res, filename, gene_data) {
    clinical = os_res
    clinical_trait <- clinical  %>%
    dplyr::select(Sample.ID,gender,vital_status,                            
                    days_to_death,days_to_last_followup,
                    person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
                    laterality,neoplasm_histologic_grade,pathologic_stage) %>%
    distinct(Sample.ID, .keep_all = TRUE)  

    #整理死亡患者的临床信息
    dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == 'Dead') %>%
    dplyr::select(-days_to_last_followup) %>%
    S4Vectors::rename(c(Sample.ID = 'Barcode',
                        gender = 'Gender',
                        vital_status = 'OS',
                        days_to_death='OS.Time',
                        person_neoplasm_cancer_status='cancer_status',
                        age_at_initial_pathologic_diagnosis = 'Age',
                        neoplasm_histologic_grade = 'Grade',
                        pathologic_stage = 'Stage')) %>%
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=OS.Time/365)

    #整理生存患者的临床信息
    alive_patient <- clinical_trait %>%
    dplyr::filter(vital_status == 'Alive') %>%
    dplyr::select(-days_to_death) %>%
    S4Vectors::rename(c(Sample.ID = 'Barcode',
                        gender = 'Gender',
                        vital_status = 'OS',
                        days_to_last_followup='OS.Time',
                        person_neoplasm_cancer_status='cancer_status',
                        age_at_initial_pathologic_diagnosis = 'Age',
                        neoplasm_histologic_grade = 'Grade',
                        pathologic_stage = 'Stage')) %>%
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=OS.Time/365)
    survival_data <- rbind(dead_patient,alive_patient)
    survival_data <- subset(survival_data,select=c(Barcode,OS ,OS.Time))
    survival_data = survival_data[, c('OS.Time', 'OS')]

    data = cbind(gene_data, survival_data)
    res <- coxph(Surv(OS.Time, OS) ~ ., data =  data)
    mul_cox <- summary(res)
    mul_HR<- round(mul_cox$coefficients[,2],2) 
    mul_Pvalue<- round(mul_cox$coefficients[,5],4) 
    mul_CI5<-round(mul_cox$conf.int[,3],2)
    mul_CI95<-round(mul_cox$conf.int[,4],2)
    mul_CI<-paste0(mul_HR,' (',mul_CI5,'-',mul_CI95,')')
    Variable<-row.names(data.frame(mul_cox$coefficients))
    mulcox_res<- data.frame(Variable,mul_HR,mul_CI5,mul_CI95,mul_CI,mul_Pvalue)
    colnames(mulcox_res)=c("Variable","HR","CI5","CI95","HR (95% CI)","Pvalue")
    write.csv(mulcox_res, paste0("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_multi_res/", filename, ".csv"))

}

setwd("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_res/")
fs = list.files('./', 'csv')
data = data[os_res$Sample.ID, ]
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    t = read.csv(x, row.names=1)
    t = t[t$Pvalue <= 0.05, ]
    print(filename)
    if (dim(t)[1] != 0) {
      gene_data = as.data.frame(data[, which(colnames(data) %in% t$Variable)])
      rownames(gene_data) = rownames(data)
      colnames(gene_data) = t$Variable
      cox_cal(os_res, filename, gene_data)
    }
})