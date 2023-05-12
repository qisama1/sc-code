library(SummarizedExperiment)
library(dplyr)
library(survival)
library(survminer)


setwd("/public/home/yuwenqi/sc-data/selected/append_ana/TCGA/PDAC/os/")
fs = list.files('./', 'os')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    marker = str_split(filename, '[_os;]',simplify = T)[,1]
    data = read.csv(x)
    data2 = read.csv(x, row.names=1)
    os_plot(data, paste0("./plot-3-7/", filename, ".pdf"), marker, data2)
})

lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    marker = str_split(filename, '[_os;]',simplify = T)[,1]
    data = read.csv(x)
    data2 = read.csv(x, row.names=1)
    os_plot2(data, paste0("./plot-5-5/", filename, ".pdf"), marker, data2)
})

lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    marker = str_split(filename, '[_os;]',simplify = T)[,1]
    data = read.csv(x)
    data2 = read.csv(x, row.names=1)
    os_plot2(data, paste0("./plot-5-5-2/", filename, ".pdf"), marker, data2)
})

os_plot = function(data, filename, marker, data2) {
    clinical = data
    clinical_trait <- clinical  %>%
    dplyr::select(X,gender,vital_status,                            
                    days_to_death,days_to_last_followup,
                    person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
                    laterality,neoplasm_histologic_grade,pathologic_stage) %>%
    distinct(X, .keep_all = TRUE)  

    #整理死亡患者的临床信息
    dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == 'Dead') %>%
    dplyr::select(-days_to_last_followup) %>%
    S4Vectors::rename(c(X = 'Barcode',
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
    S4Vectors::rename(c(X = 'Barcode',
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

    # 表达量数据
    data = data2

    marker <- marker


    data1 <- data[, which(colnames(data) %in% marker)]

    data1 <- as.data.frame(data1)
    rownames(data1) = rownames(data)
    colnames(data1) = 'gene'
    #row.names(data1) <- substr(row.names(data1),start = 1,stop = 12)
    row.names(data1) <- chartr(old='.',new = '-',x=row.names(data1))
    data1$Barcode <- row.names(data1)
    #data1 <- subset(data1,select= c(subTME1, subTME2, subTME3, subTME4,Barcode))
    dt <- merge(data1,survival_data ,by='Barcode')

    # 生存分析
    dt_OS <-  dt %>%
    dplyr::select(OS,OS.Time,gene)
    dt_OS <- na.omit(dt_OS)
    
    dt_OS[dt_OS['gene'] > quantile(dt_OS$gene, 0.7), 'type'] = 'high'
    dt_OS[dt_OS['gene'] <= quantile(dt_OS$gene, 0.3), 'type'] = 'low'


    dt_OS = subset(dt_OS, subset = type == 'high' | type == 'low')
    fit <- survfit(Surv(OS.Time,OS)~type,data=dt_OS)


    ggpar(ggsurvplot(
    fit,   data = dt_OS,   
    pval = TRUE,           
    legend.title = marker,
    xlim = c(0,5),  
    xlab = "Time in years",  
    break.time.by = 1,  
    pval.size = 8,
    legend.labs =  c("high", "Low"),   
    palette = c("#E41A1C","#377EB8")),
    font.y  = c(16, "bold"), 
    font.x  = c(16, "bold"),
    legend = "top",
    font.legend = c(16, "bold"))
    ggsave(filename = filename , device = 'pdf', width = 5, height = 5, units = 'in')

}

os_plot2 = function(data, filename, marker, data2) {
    clinical = data
    clinical_trait <- clinical  %>%
    dplyr::select(X,gender,vital_status,                            
                    days_to_death,days_to_last_followup,
                    person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
                    laterality,neoplasm_histologic_grade,pathologic_stage) %>%
    distinct(X, .keep_all = TRUE)  

    #整理死亡患者的临床信息
    dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == 'Dead') %>%
    dplyr::select(-days_to_last_followup) %>%
    S4Vectors::rename(c(X = 'Barcode',
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
    S4Vectors::rename(c(X = 'Barcode',
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

    # 表达量数据
    data = data2

    marker <- marker


    data1 <- data[, which(colnames(data) %in% marker)]

    data1 <- as.data.frame(data1)
    rownames(data1) = rownames(data)
    colnames(data1) = 'gene'
    #row.names(data1) <- substr(row.names(data1),start = 1,stop = 12)
    row.names(data1) <- chartr(old='.',new = '-',x=row.names(data1))
    data1$Barcode <- row.names(data1)
    #data1 <- subset(data1,select= c(subTME1, subTME2, subTME3, subTME4,Barcode))
    dt <- merge(data1,survival_data ,by='Barcode')

    # 生存分析
    dt_OS <-  dt %>%
    dplyr::select(OS,OS.Time,gene)
    dt_OS <- na.omit(dt_OS)
    
    dt_OS[dt_OS['gene'] > quantile(dt_OS$gene, 0.5), 'type'] = 'high'
    dt_OS[dt_OS['gene'] <= quantile(dt_OS$gene, 0.5), 'type'] = 'low'


    dt_OS = subset(dt_OS, subset = type == 'high' | type == 'low')
    fit <- survfit(Surv(OS.Time,OS)~type,data=dt_OS)


    ggpar(ggsurvplot(
    fit,   data = dt_OS,   
    pval = TRUE,           
    legend.title = marker,
    xlim = c(0,5),  
    xlab = "Time in years",  
    break.time.by = 1,  
    pval.size = 8,
    legend.labs =  c("high", "Low"),   
    palette = c("#E41A1C","#377EB8")),
    font.y  = c(16, "bold"), 
    font.x  = c(16, "bold"),
    legend = "top",
    font.legend = c(16, "bold"))
    ggsave(filename = filename , device = 'pdf', width = 5, height = 5, units = 'in')

}

os_plot = function(data, filename, marker, data2) {
    clinical = data
    clinical_trait <- clinical  %>%
    dplyr::select(X,gender,vital_status,                            
                    days_to_death,days_to_last_followup,
                    person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
                    laterality,neoplasm_histologic_grade,pathologic_stage) %>%
    distinct(X, .keep_all = TRUE)  

    #整理死亡患者的临床信息
    dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == 'Dead') %>%
    dplyr::select(-days_to_last_followup) %>%
    S4Vectors::rename(c(X = 'Barcode',
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
    S4Vectors::rename(c(X = 'Barcode',
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

    # 表达量数据
    data = data2

    marker <- marker


    data1 <- data[, which(colnames(data) %in% marker)]

    data1 <- as.data.frame(data1)
    rownames(data1) = rownames(data)
    colnames(data1) = 'gene'
    #row.names(data1) <- substr(row.names(data1),start = 1,stop = 12)
    row.names(data1) <- chartr(old='.',new = '-',x=row.names(data1))
    data1$Barcode <- row.names(data1)
    #data1 <- subset(data1,select= c(subTME1, subTME2, subTME3, subTME4,Barcode))
    dt <- merge(data1,survival_data ,by='Barcode')

    # 生存分析
    dt_OS <-  dt %>%
    dplyr::select(OS,OS.Time,gene)
    dt_OS <- na.omit(dt_OS)
    
    dt_OS[dt_OS['gene'] > quantile(dt_OS$gene, 0.7), 'type'] = 'high'
    dt_OS[dt_OS['gene'] <= quantile(dt_OS$gene, 0.3), 'type'] = 'low'


    dt_OS = subset(dt_OS, subset = type == 'high' | type == 'low')
    fit <- survfit(Surv(OS.Time,OS)~type,data=dt_OS)


    ggpar(ggsurvplot(
    fit,   data = dt_OS,   
    pval = TRUE,           
    legend.title = marker,
    xlim = c(0,5),  
    xlab = "Time in years",  
    break.time.by = 1,  
    pval.size = 8,
    legend.labs =  c("high", "Low"),   
    palette = c("#E41A1C","#377EB8")),
    font.y  = c(16, "bold"), 
    font.x  = c(16, "bold"),
    legend = "top",
    font.legend = c(16, "bold"))
    ggsave(filename = filename , device = 'pdf', width = 5, height = 5, units = 'in')

}

os_plot2 = function(data, filename, marker, data2) {
    clinical = data
    clinical_trait <- clinical  %>%
    dplyr::select(X,gender,vital_status,                            
                    days_to_death,days_to_last_followup,
                    person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
                    laterality,neoplasm_histologic_grade,pathologic_stage) %>%
    distinct(X, .keep_all = TRUE)  

    #整理死亡患者的临床信息
    dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == 'Dead') %>%
    dplyr::select(-days_to_last_followup) %>%
    S4Vectors::rename(c(X = 'Barcode',
                        gender = 'Gender',
                        vital_status = 'OS',
                        days_to_death='OS.Time',
                        person_neoplasm_cancer_status='cancer_status',
                        age_at_initial_pathologic_diagnosis = 'Age',
                        neoplasm_histologic_grade = 'Grade',
                        pathologic_stage = 'Stage')) %>%
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=OS.Time/12)

    #整理生存患者的临床信息
    alive_patient <- clinical_trait %>%
    dplyr::filter(vital_status == 'Alive') %>%
    dplyr::select(-days_to_death) %>%
    S4Vectors::rename(c(X = 'Barcode',
                        gender = 'Gender',
                        vital_status = 'OS',
                        days_to_last_followup='OS.Time',
                        person_neoplasm_cancer_status='cancer_status',
                        age_at_initial_pathologic_diagnosis = 'Age',
                        neoplasm_histologic_grade = 'Grade',
                        pathologic_stage = 'Stage')) %>%
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=OS.Time/12)
    survival_data <- rbind(dead_patient,alive_patient)
    survival_data <- subset(survival_data,select=c(Barcode,OS ,OS.Time))

    # 表达量数据
    data = data2

    marker <- marker


    data1 <- data[, which(colnames(data) %in% marker)]

    data1 <- as.data.frame(data1)
    rownames(data1) = rownames(data)
    colnames(data1) = 'gene'
    #row.names(data1) <- substr(row.names(data1),start = 1,stop = 12)
    row.names(data1) <- chartr(old='.',new = '-',x=row.names(data1))
    data1$Barcode <- row.names(data1)
    #data1 <- subset(data1,select= c(subTME1, subTME2, subTME3, subTME4,Barcode))
    dt <- merge(data1,survival_data ,by='Barcode')

    # 生存分析
    dt_OS <-  dt %>%
    dplyr::select(OS,OS.Time,gene)
    dt_OS <- na.omit(dt_OS)
    
    dt_OS[dt_OS['gene'] > quantile(dt_OS$gene, 0.5), 'type'] = 'high'
    dt_OS[dt_OS['gene'] <= quantile(dt_OS$gene, 0.5), 'type'] = 'low'


    dt_OS = subset(dt_OS, subset = type == 'high' | type == 'low')
    fit <- survfit(Surv(OS.Time,OS)~type,data=dt_OS)


    ggpar(ggsurvplot(
    fit,   data = dt_OS,   
    pval = TRUE,           
    legend.title = marker,
    #xlim = c(0,5),  
    xlab = "Time in months",  
    break.time.by = 1,  
    pval.size = 8,
    legend.labs =  c("high", "Low"),   
    palette = c("#E41A1C","#377EB8")),
    font.y  = c(16, "bold"), 
    font.x  = c(16, "bold"),
    legend = "top",
    font.legend = c(16, "bold"))
    ggsave(filename = filename , device = 'pdf', width = 5, height = 5, units = 'in')

}