library(qs)
library(tidyverse)
data = qread("/public/home/yuwenqi/targers_work/data_concat/data_used_protein.qs")
sample = read.csv("/public/home/yuwenqi/targers_work/data_concat/concated_sample.csv", row.names = 1)
data[is.na(data)] = 0   
data['cancer_type'] = sample[rownames(data), 'cancer_type']
data['cancer'] = sample[rownames(data), 'Tissue']
cancers = read.csv("/public/home/yuwenqi/targers_work/data_concat/check/cancers.csv")

setwd("/public/home/yuwenqi/works/zhangzimei/filter/")

library(stringr)

protein_genes = read.csv("/public/home/yuwenqi/targers_work/data_concat/val_filter/id_name_protein_coding22.txt", sep = ',')

for (cancer in cancers$Tissue) {
    path = paste0("/public/home/yuwenqi/targers_work/data_concat/val_filter/normal_filter/", cancer, "/")
    setwd(path)
    fs=list.files('./','csv')
    cancer_data = subset(data, cancer == cancer)
    tumor_data = subset(cancer_data, cancer_type == 'Tumor')
    normal_data = subset(cancer_data, cancer_type == 'Normal')
    for (num in fs) {
        num_str = str_split(num,'.csv',simplify = T)[, 1]
        num_cur = as.numeric(num_str)
        gene_selected = data.frame()
        genes = read.csv(num)
        for (gene in genes$gene[genes$gene %in% colnames(data)]){
            cnt = 0
            if (gene %in% colnames(data)) {
                tumor_used = tumor_data[, c(gene, 'cancer_type')]
                normal_used = normal_data[, c(gene, 'cancer_type')]
                colnames(tumor_used) = c('gene', 'cancer_type')
                colnames(normal_used) = c('gene', 'cancer_type')
                if (!is.na(wilcox.test(tumor_used$gene, normal_used$gene)$p.value) && wilcox.test(tumor_used$gene, normal_used$gene)$p.value < 0.05 && mean(tumor_used[, 'gene']) > 2 * mean(normal_used[, 'gene'])) {
                    # for (cancer in cancers) {
                    #     cancer_data = subset(data, cancer == cancer)
                    #     tumor_data = subset(cancer_data, cancer_type == 'Tumor')
                    #     normal_data = subset(cancer_data, cancer_type == 'Normal')
                    #     tumor_used = tumor_data[, c(gene, 'cancer_type')]
                    #     normal_used = normal_data[, c(gene, 'cancer_type')]
                    #     colnames(tumor_used) = c('gene', 'cancer_type')
                    #     colnames(normal_used) = c('gene', 'cancer_type')
                    #     if (!is.na(wilcox.test(tumor_used$gene, normal_used$gene)$p.value) && wilcox.test(tumor_used$gene, normal_used$gene)$p.value < 0.05 && mean(tumor_used[, 'gene']) > 2 * mean(normal_used[, 'gene'])) {
                    #         cnt = cnt + 1
                    #     }
                    # }
                        t = data.frame(gene)
                        rownames(t) = c(gene)
                        colnames(t) = c('genes')
                        gene_selected = rbind(gene_selected, t)
                        cnt = 1
                }
                print(paste0(gene, " ", cnt))
            }
            # if (cnt >= num_cur) {
            #     t = data.frame(gene)
            #     rownames(t) = c(gene)
            #     colnames(t) = c('genes')
            #     gene_selected = rbind(gene_selected, t)
            # }
            
        }
        
        save_path = paste0("/public/home/yuwenqi/targers_work/data_concat/val_filter/tumor_normal/normal/", cancer, "/")
        # save_path = "F://Rstdio/work"
        # num = "14.csv"
        # save_path + num =  "F://Rstdio/work/14.csv"
        dir.create(save_path)
        write.csv(gene_selected, paste0(save_path, num))
    }
}

