library(CellChat)
library(patchwork)
library(Seurat)
library(qs)
options(stringsAsFactors = FALSE)

scRNA = qread("/public/home/yuwenqi/sc-data/selected/9/all.qs")
meta = read.csv("/public/home/yuwenqi/sc-data/selected/9/all3.csv", row.names=1)
scRNA = scRNA[, rownames(meta)]
scRNA[['ident']] = meta['ident']

module = read.csv("/public/home/yuwenqi/sc-data/selected/9/module/module3.csv", encoding = 'GBK')
scRNA_m1 = scRNA[, which(scRNA$ident %in% module$module3)]

data.input = scRNA_m1@assays$RNA@data
meta.data <- subset(scRNA_m1@meta.data, select = c("orig.ident", "ident"))

meta.data$sub_clusters = factor(meta.data$ident)

### 1.2 Create a CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "sub_clusters")
                           
### 1.4 加载CellChat受配体数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat,population.size = F)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
head(df.net)
write.csv(df.net, "/public/home/yuwenqi/sc-data/selected/9/module/cellchat/module3_df_net.csv")
#以通路为单位提取通讯信息
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
qsave(cellchat, "/public/home/yuwenqi/sc-data/selected/CRC/workspace/CCI/cellchat/cellchat.qs")