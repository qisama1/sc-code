library(readr)
pathway <- read_delim("/public/home/yuwenqi/sc-data/selected/new-RNAseq/pathway.csv", ",", 
                      escape_double = FALSE, trim_ws = TRUE)
                    
pathway <- as.data.frame(pathway)
pathway_list <- lapply(pathway, function(x) {
  unique(na.omit(x)) 
})

tpm = read.csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/all-tpm_unstranded.csv", row.names=1)
brca_meta = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/BRCA_sample_meta2.csv", row.names=1)
pdac_meta = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/pdac_sample_meta2.csv", row.names=1)
esca_meta = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ESCA_sample_meta2.csv", row.names=1)
crc_meta = read.csv("/public/home/yuwenqi/sc-data/selected/new-RNAseq/COAD/coad_meta2.csv", row.names=1)
id_name = read.csv("/public2022/chenzixi/research/TCGA_data/new-RNAseq/id_name.csv", row.names=1)
tpm = t(tpm)
colnames(tpm) = id_name[colnames(tpm), 'Gene_name']
tpm = t(tpm)

brca_tpm = tpm[, colnames(brca_meta)]
pdac_tpm = tpm[, colnames(pdac_meta)]
esca_tpm = tpm[, colnames(esca_meta)]
crc_tpm = tpm[, colnames(crc_meta)]

brca_ssgsea_score = gsva(brca_tpm, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(brca_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/ssgsea/brca_ssgsea.csv")

pdac_ssgsea_score = gsva(pdac_tpm, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(pdac_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/ssgsea/pdac_ssgsea.csv")

esca_ssgsea_score = gsva(esca_tpm, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(esca_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ssgsea/esca_ssgsea.csv")

crc_ssgsea_score = gsva(crc_tpm, pathway_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
write.csv(crc_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/COAD/gsva/crc_ssgsea.csv")


es <- gsva(as.matrix(exp), geneset, kcdf='Gaussian', method = 'gsva', parallel.sz=4, verbose = TRUE)

brca_ssgsva_score = gsva(brca_tpm, pathway_list, kcdf='Gaussian', method = 'gsva', verbose = TRUE)   # signature 'matrix,list'
write.csv(brca_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/BRCA/ssgsea/brca_ssgsva.csv")

pdac_ssgsva_score = gsva(pdac_tpm, pathway_list, kcdf='Gaussian', method = 'gsva', verbose = TRUE)   # signature 'matrix,list'
write.csv(pdac_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/PDAC/ssgsea/pdac_ssgsva.csv")

esca_ssgsva_score = gsva(esca_tpm, pathway_list, kcdf='Gaussian', method = 'gsva', verbose = TRUE)   # signature 'matrix,list'
write.csv(esca_ssgsea_score, "/public/home/yuwenqi/sc-data/selected/new-RNAseq/ESCA/ssgsea/esca_ssgsva.csv")
