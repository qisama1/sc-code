## load data
library(Seurat)
library(dplyr)
library(monocle3)
options(stringsAsFactors=FALSE)
library(reticulate)

glm.predict <-
  function(train.data, train.group, downsample = FALSE, sample.cells = 0, genes.used = NA, test.data, test.group, alpha = 0.99, nfolds = 10) {
    ## Calculate the similarities of the train data and test data.
    ##
    ## Args:
    #' @train.data: A train data matrix with each cell in a column and each gene
    #' in a row.
    #' @train.group: A vector with the same length as the column of train.data.
    #' @downsample: Whether to sample cells in each cluster to the minimum cluster size.
    #' @sample.cells: Sample cells in each group of cells in train data, if 0 do not
    #' sample cells.
    #' @genes.used: Use a subset of genes in both the train and test data.
    #' @test.data: A test data matrix with each cell in a column and each gene
    #' in a row.
    #' @test.group: A vector with the same length as the column of train.data.
    #' @alpha: The elasticnet mixing parameter, with 0≤α≤1, passed to cv.glmnet.
    #' @nfolds: Number of folds, passed to cv.glmnet.
    ##
    ## Returns:
    ## The probability of each cell in the test.data to be predicted as each group.
    require(glmnet)
    require(ComplexHeatmap)
    glm.fits <- list()
    glm.predict <- list()
    if (length(genes.used) > 1) {
      train.data <- train.data[genes.used,]
      test.data <- test.data[genes.used,]
      if (length(genes.used) <= 50) {
        cat("There were less than 50 features used in the training data!\n")
      }
    }
    if (sample.cells == 0 & downsample) {
      sample.cells <- max(50, min(table(train.group)))
    }
    if (sample.cells > 0) {
      ngroup <- length(unique(train.group))
      if (ncol(train.data) >= sample.cells * ngroup) {
        cells_used <- c()
        for (groupi in sort(unique(train.group))) {
          if (length(which(train.group == groupi)) > sample.cells) {
            cells_used <-
              c(cells_used, sample(which(train.group == groupi), sample.cells))
          } else{
            cells_used <- c(cells_used, which(train.group == groupi))
          }
        }
        train.data <- train.data[, cells_used]
        train.group <- train.group[cells_used]
      }
    }
    for (groupi in sort(unique(train.group))) {
      fac <-  factor(train.group == groupi)
      glm.fits[[groupi]] <-
        cv.glmnet(x = t(train.data), fac, offset = getPopulationOffset(fac), 
                  family = 'binomial', intercept = FALSE, 
                  alpha = alpha, nfolds = nfolds, type.measure = 'class'
        )
      glm.predict[[groupi]] <-
        predict(
          object = glm.fits[[groupi]],
          newx = t(test.data),
          newoffset = rep(0, ncol(test.data)),
          s = 'lambda.min'
        )
    }
    glm.predict.df <- data.frame(do.call(cbind, glm.predict))
    colnames(glm.predict.df) <- sort(unique(train.group))
    glm.predict.df.prob <- (1 + exp(-glm.predict.df)) ** -1
    glm.cluster <-
      colnames(glm.predict.df.prob)[apply(glm.predict.df.prob, 1, which.max)]
    glm.predict.mean <-
      apply(glm.predict.df, 2, function(e)
        sapply(split(e, test.group), mean))
    glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
    heatmap <- Heatmap(
      t(glm.predict.mean.prob),
      name = 'Predicted\nSimilarity',
      column_title = 'test data',
      row_title = 'train data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
    return(
      list(
        test.group = test.group,
        logits = glm.predict.df,
        probability = glm.predict.df.prob,
        cluster = glm.cluster,
        heatmap = heatmap
      )
    )
  }

getPopulationOffset = function(y) {
  ## Calculate the offset value used in glm.predict.
  if (!is.factor(y))
    y = factor(y)
  if (length(levels(y)) != 2)
    stop("y must be a two-level factor")
  off = sum(y == levels(y)[2]) / length(y)
  off = log(off / (1 - off))
  return(rep(off, length(y)))
}

data1 = qread("/public/home/yuwenqi/sc-data/selected/10/end.qs")  #dataset10
data2 = qread("/public/home/yuwenqi/sc-data/selected/9/endo.qs") #dataset9

meta1 = read.csv("/public/home/yuwenqi/sc-data/selected/10/end_info.csv", row.names=1)
meta2 = read.csv("/public/home/yuwenqi/sc-data/selected/9/end_info.csv", row.names=1)

data1$ident = meta1['ident']
data2$ident = meta2['ident']

ada1_df = data1$ident
ada2_df = data2$ident

gene1 <- VariableFeatures(data1)
gene2 <- VariableFeatures(data2)


ada1_exp <- as.matrix(data1@assays$RNA@data)
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- as.matrix(data2@assays$RNA@data)
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2))


res <- glm.predict(ada2_exp, ada2_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada1_exp, ada1_df, alpha = 0.99, nfolds = 10)

library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("#e9e9e9","white", "red"))
glm.predict.mean <-
      apply(res$logits, 2, function(e)
        sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
plot <- Heatmap(
      t(glm.predict.mean.prob),
	  col = col_fun,
      name = 'Predicted\nSimilarity',
      column_title = 'data10',
      row_title = 'data9',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
pdf("/public/home/yuwenqi/sc-data/selected/9-10-combined/combine_end.pdf")
plot
dev.off()