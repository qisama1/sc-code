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

coef <- summary(coxm)[["coefficients"]]
data_risk <- as.data.frame((tcga_data))
data_risk <- subset(data_risk,select=gene_list)
k <- length(gene_list)
i <- 1
riskScore <- 0
while (i <= k) {
  riskScore <- riskScore+coef[i,1]*data_risk[,i]
  i <- i + 1
}
data_risk$riskScore <- riskScore

risk_data = as.data.frame(data_risk$riskScore)
colnames(risk_data) = c('risk')
rownames(risk_data) = rownames(data_risk)
res_path = paste0("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/cluster-an/cox/module_cluster_genes/subTME1/risk2/", "test.csv")

write.csv(risk_data, res_path)
