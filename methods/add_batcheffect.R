library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

add_batch_effect <- function(embryo_sce, batch_col="synth_samples", norm_sd=0.5){
  cellids_sample <- split(embryo_sce$cell, embryo_sce[[batch_col]])
  X_pca <- reducedDim(embryo_sce, "HARMONY")
  X_pca_batch <- X_pca

  for (b in names(cellids_sample)){
    batch_effect <- rnorm(ncol(X_pca), mean=0, sd = norm_sd)
    X_pca_batch[cellids_sample[[b]],] <- t(apply(X_pca_batch[cellids_sample[[b]],], 1, function(x) x + batch_effect))
  }
  
  reducedDim(embryo_sce, "pca_batch") <- X_pca_batch
  embryo_sce
}

# main
data <- readRDS("/.../ATACcydar_test.rds")

coldata <- read.csv("/.../CRClabel.csv")

colData(data)$cell <- data.frame(coldata)$rowname
colData(data)$synth_samples <- data.frame(coldata)$synth_samples

cond_probability<-coldata[,c("Condition1_prob","Condition2_prob")]
colnames(cond_probability) = paste0("Condition", 1:2)

synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
colData(data)$synth_labels <-synth_labels

#sd=0.75, 1, 1.25, 1.5
data_be <- add_batch_effect(data, batch_col="synth_labels", norm_sd=1.5)



