data = read.csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/coad_tcga.csv", row.names = 1)
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/COAD/os/os_meta.csv")

setwd("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_multi_res/")
fs = list.files('./', 'csv')
lapply(fs, function(x) {
    print(x)
    filename = str_split(x,'[.;]',simplify = T)[,1]
    cox_path = paste0("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_res/", x)
    cox = read.csv(cox_path, row.names=1)
    t = read.csv(x, row.names=1)
    t = t[t$Pvalue <= 0.05, ]
    variables = as.data.frame(cox[which(cox$Variable %in% t$Variable), ])
    if (dim(variables)[1] != 0) {
        cur = variables[1, ]
        data$risk = data[cur$Variable] * cur$HR
        for (i in 1:nrow(variables)) {
            cur = variables[i, ]
            data$risk = data$risk + data[cur$Variable] * cur$HR
        }
        data = data[os_res$Sample.ID, ]
        res_path = paste0("/public/home/yuwenqi/sc-data/selected/CRC/workspace/cluster-an/cox/crc/module_cluster_genes/subTME1/cox_risk/" , x)
        risk_data = data$risk
        colnames(risk_data) = c('risk')
        write.csv(risk_data, res_path)
    }
    
})