library("survival")
#以R包survival中的肺癌数据为例
data = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/coad_tcga.csv", row.names = 1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/COAD/os/os_meta.csv")

tcga_data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/tcga_tpm_log.csv", row.names=1)
tcga_meta = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/BRCA_sample_meta2.csv", row.names=1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/os/os_res.csv", row.names=1)
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

    data = cbind(gene_data, survival_data)
    y<- Surv(data$OS.Time, data$OS)
    Unicox_model<- function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=data)#
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],2)#指定小数点后数位
    Pvalue <- round(cox1$coefficients[,5],3)
    CI5 <-round(cox1$conf.int[,3],2)
    CI95 <-round(cox1$conf.int[,4],2)
    CI<- paste0(HR," (",CI5,"-",CI95,")")
    Variable=row.names(data.frame(cox1$coefficients))
    Unicox_model<- data.frame('Variable' = Variable,
                                'HR' = HR,
                                'HR(95% CI)' = CI,
                                'Pvalue' = Pvalue)
    return(Unicox_model)}  
    variable<- colnames(gene_data)#设置变量
    Unicox <- lapply(variable, Unicox_model)#循环输入
    Unicox<- ldply(Unicox,data.frame)#将list转变为数据框
    colnames(Unicox)=c('Variable', 'HR', 'HR (95% CI)','Pvalue')
    write.csv(Unicox, paste0("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_res/", filename, ".csv"))
}

setwd("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/")
fs = list.files('./', 'csv')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    genes = read.csv(x, row.names=1)
    gene_data = data[, which(colnames(data) %in% genes$gene)]
    cox_cal(os_res, filename, gene_data[os_res$Sample.ID, ])
})


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
        

    data = cbind(gene_data, survival_data)
    y<- Surv(data$OS.Time, data$OS)
    Unicox_model<- function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=data)#
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],2)#指定小数点后数位
    Pvalue <- round(cox1$coefficients[,5],3)
    CI5 <-round(cox1$conf.int[,3],2)
    CI95 <-round(cox1$conf.int[,4],2)
    CI<- paste0(HR," (",CI5,"-",CI95,")")
    Variable=row.names(data.frame(cox1$coefficients))
    Unicox_model<- data.frame('Variable' = Variable,
                                'HR' = HR,
                                'HR(95% CI)' = CI,
                                'Pvalue' = Pvalue)
    return(Unicox_model)}  
    variable<- colnames(gene_data)#设置变量
    Unicox <- lapply(variable, Unicox_model)#循环输入
    Unicox<- ldply(Unicox,data.frame)#将list转变为数据框
    colnames(Unicox)=c('Variable', 'HR', 'HR (95% CI)','Pvalue')
    write.csv(Unicox, paste0("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_res/", filename, ".csv"))
}
