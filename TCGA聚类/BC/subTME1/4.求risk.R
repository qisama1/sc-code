tcga_data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/split_data/BRCA_data.csv", row.names=1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/os/os_res.csv")

os_res=os_res[os_res$vital_status %in% c('Alive','Dead'),]
os_res$days_to_last_follow_up[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
os_res$os = ifelse(os_res$vital_status =='Dead', 1, 0)
os_res = os_res[os_res$days > 0,]
tcga_data = tcga_data[os_res$Sample.ID, ]

setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME1/multi_cox_uni_res/")
fs = list.files('./', 'csv')
lapply(fs, function(x) {
    print(x)
    filename = str_split(x,'[.;]',simplify = T)[,1]
    cox_path = paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME1/cox_res/", x)
    cox = read.csv(cox_path, row.names=1)
    t = read.csv(x, row.names=1)
    t = t[t$Pvalue <= 0.05, ]
    variables = as.data.frame(cox[which(cox$Variable %in% t$Variable), ])
    if (dim(variables)[1] != 0) {
        cur = variables[1, ]
        tcga_data$risk = tcga_data[cur$Variable] * cur$HR
        for (i in 1:nrow(variables)) {
            cur = variables[i, ]
            tcga_data$risk = tcga_data$risk + tcga_data[cur$Variable] * cur$HR
        }
        tcga_data = tcga_data[os_res$Sample.ID, ]
        res_path = paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME1/risk/" , x)
        risk_data = tcga_data$risk
        colnames(risk_data) = c('risk')
        write.csv(risk_data, res_path)
    }
    
})