library(SummarizedExperiment)
library(dplyr)
library(survival)
library(survminer)

vital_status， day to death

clinical = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ESCA_meta.csv")
clinical_trait <- clinical  %>%
  dplyr::select(bcr_patient_barcode,gender,vital_status,                            
                days_to_death,days_to_last_followup,
                person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
                laterality,neoplasm_histologic_grade,pathologic_stage) %>%
  distinct( bcr_patient_barcode, .keep_all = TRUE)  

#整理死亡患者的临床信息
dead_patient <- clinical_trait  %>%
  dplyr::filter(vital_status == 'Dead') %>%
  dplyr::select(-days_to_last_followup) %>%
  S4Vectors::rename(c(bcr_patient_barcode = 'Barcode',
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
  rename(c(bcr_patient_barcode = 'Barcode',
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
data = clinical

marker <- c('module1', 'module2', 'module3')


data1 <- data[,which(colnames(data) %in% marker)]

data1 <- as.data.frame(data1)
#colnames(data1) = 'ERBB2'
row.names(data1) <- substr(row.names(data1),start = 1,stop = 12)
row.names(data1) <- chartr(old='.',new = '-',x=row.names(data1))
data1$Barcode <- row.names(data1)
data1 <- subset(data1,select= c(ERBB2,Barcode))
dt <- merge(data1,survival_data ,by='Barcode')

# 生存分析
dt_OS <-  dt %>%
  dplyr::select(OS,OS.Time,ERBB2)
dt_OS <- na.omit(dt_OS)
dt_OS['ERBB2'] <- ifelse(dt_OS['ERBB2'] > median(dt_OS[,'ERBB2']), 'high', 'low')


fit <- survfit(Surv(OS.Time,OS)~ERBB2,data=dt_OS)

ggpar(ggsurvplot(
  fit,   data = dt_OS,   
  pval = TRUE,           
  legend.title = "ERBB2",
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
ggsave(filename = "/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/survival/ERBB2_os.pdf", device = 'pdf', width = 5, height = 5, units = 'in')

exp = t(as.matrix(as.data.frame(scRNA['ERBB2']@assays$RNA@data)))
exp <- na.omit(exp)
scRNA[['ERBB2_status']] <- ifelse(exp[,'ERBB2'] > median(exp[,'ERBB2']), 'high', 'low')

ifelse(exp['ERBB2'] > median(exp[, 'ERBB2']), 'high', 'low')