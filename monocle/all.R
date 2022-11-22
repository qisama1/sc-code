library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)

setwd("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/monocle/all/")


## Runing in r-4.1.3
meta = read.csv("/public/home/yuwenqi/sc-data/selected/Breast/workspace2/BC3-CCI/bc3_meta_all.csv", row.names=1)
pbmc[['cell_type']] = meta['sub_clusters']