library(tidyverse)
library(survival)
library(survminer)
library(glmnet)
library(rms)
library(survivalROC)
library(timeROC)
library(forestplot)
library(plyr)

tcga_data = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/new_process/split_data/BRCA_data.csv", row.names=1)
tcga_meta = t(read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/BRCA_sample_meta2.csv", row.names=1))
os_res = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/os/os_res.csv")

os_res=os_res[os_res$vital_status %in% c('Alive','Dead'),]
os_res$days_to_last_follow_up[is.na(os_res$days_to_last_followup)] = 0 #is.na()用于返回是否为缺失值
os_res$days_to_death[is.na(os_res$days_to_death)] = 0   
os_res$days<-ifelse(os_res$vital_status=='Alive',os_res$days_to_last_followup,os_res$days_to_death)
os_res$os = ifelse(os_res$vital_status =='Dead', 1, 0)
os_res = os_res[os_res$days > 0,]

lasso_cal = function(os_res, filename, unicox) {
    unicox = unicox[unicox$Pvalue <= 0.05, ]
    if (dim(unicox[1]) > 1) {
        X <- as.matrix(tcga_data[os_res$Sample.ID, which(colnames(tcga_data) %in% unicox$Variable)])
        Y <- data.matrix(Surv(os_res$days,os_res$os))
        f1 = glmnet(X, Y, family="cox", nlambda=100, alpha=1)
        
        pdf(paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/lasso/" , filename, "_" ,"lambda_plot_f1.pdf"),width=15,height=15)
        plot(f1, xvar="lambda",label=TRUE)
        dev.off()

        pdf(paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/lasso/" , filename, "_" ,"lambda_plot_cv.pdf"),width=15,height = 15)
        fitcv <- cv.glmnet(X,Y,family="cox",alpha=1,nfolds=10)
        plot(fitcv)
        dev.off()

        coef(fitcv, s="lambda.min")
        gene <- coef(fitcv, s="lambda.min")
        Active.Index <- which(as.numeric(gene)!= 0)
        active.coefficients <- as.numeric(gene)[Active.Index]
        genes = as.data.frame(gene_list <- rownames(gene)[Active.Index])
        colnames(genes) = c('gene')
        write.csv(genes, paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/lasso/" , filename, "_" ,"genes.csv"))

        X <- as.data.frame(X)
        Y <- as.data.frame(Y)
        X <- subset(X,select = gene_list)
        X$time <- Y$time
        X$status <- Y$status
        coxm <- coxph(Surv(time,status)~.,data=X)
        sum_muticox <- summary(coxm)
        coef <- sum_muticox[["coefficients"]]
        p_max <- max(coef[,5])
        
        while(p_max > 0.05){
            coef <- as.data.frame(coef)
            colnames(coef)[5] <- "P.value"
            coef$name <- row.names(coef)
            coef <- arrange(coef,desc(P.value))[-1,]
            formula_for_multicox <- as.formula(paste0('Surv(time,status)~', paste(coef$name, sep = "",collapse = '+')))
            coxm <- coxph(formula_for_multicox,data=X)
            sum_muticox <- summary(coxm)
            coef <- sum_muticox[["coefficients"]]
            p_max <- max(coef[,5])
        }
        gene_list <- row.names(coef)
        muticox <- summary(coxm)[["coefficients"]]
        res_path = paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/risk2/", filename, "_multicox.csv")
        write.csv(muticox, res_path)
        coef <- summary(coxm)[["coefficients"]]
        data_risk <- as.data.frame((tcga_data))
        data_risk <- subset(data_risk,select=gene_list)
        k <- length(gene_list)
        if (k != 0) {
            i <- 1 
            riskScore <- 0
            while (i <= k) {
            riskScore <- riskScore+coef[i,2]*data_risk[,i]
            i <- i + 1
            }
            data_risk$riskScore <- riskScore

            risk_data = as.data.frame(data_risk$riskScore)
            colnames(risk_data) = c('risk')
            rownames(risk_data) = rownames(data_risk)
            res_path = paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/risk2/", filename, "_risk.csv")

            write.csv(risk_data, res_path)
        }
            
    }
    
}


setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME2/cox_res/")

fs = list.files('./', 'csv')
lapply(fs, function(x) {
    filename = str_split(x,'[.;]',simplify = T)[,1]
    unicox = read.csv(x)
    # gene_data = data[, which(colnames(data) %in% genes$gene)]
    lasso_cal(os_res, filename, unicox)
})