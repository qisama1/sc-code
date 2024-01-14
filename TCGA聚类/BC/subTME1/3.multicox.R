tcga_data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/split_data/BRCA_data.csv", row.names=1)
tcga_meta = t(read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/BRCA_sample_meta2.csv", row.names=1))
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/os/os_res.csv")

os_res=os_res[os_res$vital_status %in% c('Alive','Dead'),]
os_res$days_to_last_follow_up[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
os_res$os = ifelse(os_res$vital_status =='Dead', 1, 0)
os_res = os_res[os_res$days > 0,]

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
    write.csv(mulcox_res, paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME1/multi_cox_res/", filename, ".csv"))

}


setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME1/lasso/")
fs = list.files('./', 'genes.csv')
tcga_data = tcga_data[os_res$Sample.ID, ]
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    t = read.csv(x, row.names=1)
    print(filename)
    if (dim(t)[1] != 0) {
      gene_data = as.data.frame(tcga_data[, which(colnames(tcga_data) %in% t$gene)])
      rownames(gene_data) = rownames(tcga_data)
      colnames(gene_data) = t$gene
      cox_cal(os_res, filename, gene_data)
    }
})