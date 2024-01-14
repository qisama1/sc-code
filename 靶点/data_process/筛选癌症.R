cancers = c('ESCA', 'UCEC', 'HNSC', 'KIRP', 'OV',
'LUSC', 'STAD', 'UCS', 'THCA', 'MESO', 'CESC', 'PRAD',
'COAD', 'SARC', 'LIHC', 'BRCA', 'READ',
'KICH', 'LUAD', 'BLCA', 'PAAD', 'CHOL', 'UVM', 'THYM',
'TGCT', 'ACC', 'GBM', 'SKCM', 'KIRC')

gene_list = list()
all_cancers = data.frame()
for (c in cancers) {
    cur_data = subset(data, cancer == c)
    print(table(cur_data$cancer_type))
    if (table(cur_data$cancer_type)['Normal'] > 20) {
        t = data.frame(c)
        rownames(t) = c(c)
        colnames(t) = c('cancers')
        all_cancers = rbind(all_cancers, t)
    }   
}

all_cancers = data.frame()
for (c in cancers) {
    cur_data = subset(sample, Tissue == c)
    print(table(cur_data$cancer_type))
    if (length(table(cur_data$cancer_type)) == 2) {
        if (table(cur_data$cancer_type)['Normal'] > 20) {
            t = data.frame(c)
            rownames(t) = c(c)
            colnames(t) = c('cancers')
            all_cancers = rbind(all_cancers, t)
        }   
    }
    
}

